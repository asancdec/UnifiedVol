// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#include "Math/Interpolation.hpp"
#include "Core/Matrix/Functions.hpp"
#include "Core/Functions.hpp"
#include "Utils/Aux/Errors.hpp"

#include <array>       
#include <cmath>      
#include <cstddef>     
#include <limits>    
#include <utility>     

namespace uv::models::svi
{
    template <::nlopt::algorithm Algo>
    Vector<Params> calibrate(const Vector<Real>& tenors,
        const core::Matrix<Real>& kMatrix,
        const core::Matrix<Real>& wMMatrix,
        const opt::nlopt::Optimizer<4, Algo>& prototype)
    {
        // ---------- Validate inputs ----------

        detail::validateInputs(tenors, kMatrix, wMMatrix);

        // ---------- Extract data ----------
        
        // NLopt API accepts double type only
        const core::Matrix<double> kMatrixD{ core::convertMatrix<double>(kMatrix)};
        const core::Matrix<double> wMMatrixD{ core::convertMatrix<double>(wMMatrix)};

        const std::size_t numTenors{ tenors.size()};
        const std::size_t numStrikes{ kMatrixD.cols()};

        // ---------- Initialize params container ----------

        Vector<Params> params;
        params.reserve(numTenors);

        // ---------- Calibrate per slice ----------

        for (std::size_t i = 0; i < numTenors; ++i)
        {
            const Params* prev = (i == 0) ? nullptr : &params.back();

            params.emplace_back
            (
                detail::calibrateSlice(
                    tenors[i],
                    kMatrixD[i],
                    wMMatrixD[i],
                    prototype,
                    prev,
                    numStrikes
                )
            );
        }

        return params;
    }
} // namespace uv::models::svi

namespace uv::models::svi::detail
{
    struct ObjectiveContexts
    {
        const double* k;      // Pointer to log-forward moneyness
        const double* wM;     // Pointer to market total variance
        std::size_t   n;      // Number of points
        const double  atmWK;  // Atm variance
    };

    struct CalendarContexts
    {
        double k;       // Log-forward moneyness
        double prevWk;  // Total variance of the previous slice
        double eps;     // Epsilon value
        double atmWK;   // Atm variance
    };

    struct ConvexityContexts
    {
        double k;       // Log-forward moneyness 
        double atmWK;   // Atm variance
    };

    struct GkCache
    {
        double x;             // x := k-m
        double R;             // R:= sqrt(x^2 + sigma^2)
        double invR;          // invR := 1 / R
        double wk;            // w(k) = a + b*(rho*x + R)
        double wkD1;          // w'(k) = b * (rho + x/R)
        double wkD1Squared;   // w'(k)^2
        double invRCubed;     // 1/(R^3)
        double invR5;         // 1/(R^5)
        double sigmaSquared;  // sigma^2
        double wkD2;          // w''(k) = b * sigma^2 / R^3
        double A;             // A := 1 - k * w'/(2 * w)
        double B;             // B := 1/w(k) + 1/4
        double wkInv;         // 1/w(k)
        double wkSquaredInv;  // 1/w(k)^2
        double R0;            // R0 := sqrt(m^2 + sigma^2)
        double invR0;         // 1/R0

        explicit GkCache(double a, double b, double rho, double m, double sigma, double k) noexcept
        {
            R0 = std::hypot(m, sigma);
            invR0 = 1.0 / R0;
            x = k - m;
            R = std::hypot(x, sigma);
            invR = 1.0 / R;
            wk = std::fma(b, (rho * x + R), a);
            wkInv = 1.0 / wk;
            wkSquaredInv = wkInv * wkInv;
            wkD1 = b * (rho + x * invR);
            wkD1Squared = wkD1 * wkD1;
            invRCubed = invR * invR * invR;
            invR5 = invRCubed * invR * invR;
            sigmaSquared = sigma * sigma;
            wkD2 = b * sigmaSquared * invRCubed;
            A = 1.0 - 0.5 * k * wkD1 * wkInv;
            B = wkInv + 0.25;
        }
    };

    template <::nlopt::algorithm Algo>
    Params calibrateSlice(
        const Real T,
        std::span<const double> kSlice,
        std::span<const double> wKSlice,
        const opt::nlopt::Optimizer<4, Algo>& prototype,
        const Params* prevParams,         
        const std::size_t numStrikes
    )
    {
        // ---------- Calculate  ----------

        const double atmWK{ math::interp::pchipInterp<double>(0.0, kSlice, wKSlice) };
        const double logKFMin{ core::minValue(kSlice) };
        const double logKFMax{ core::maxValue(kSlice) };

        // ---------- Fresh optimizer  ----------

        opt::nlopt::Optimizer optimizer{ prototype.fresh() };
        optimizer.setUserValue(atmWK);

        // ---------- Bounds and init guess ----------

        optimizer.setGuessBounds(
            initGuess(),
            lowerBounds(logKFMin),
            upperBounds(logKFMax)
        );

        // ---------- Standard constraints ----------

        addWMinConstraint(optimizer);
        addMinSlopeConstraint(optimizer);
        addMaxSlopeConstraint(optimizer);

        // ---------- Calendar constraints ----------

        Vector<CalendarContexts> calCtx;

        if (prevParams)
        {
            calCtx.resize(numStrikes);

            // Extract data
            const Params& prev{ *prevParams };
            const double a0{ double(prev.a) };
            const double b0{ double(prev.b) };
            const double rho0{ double(prev.rho) };
            const double m0{ double(prev.m) };
            const double sigma0{ double(prev.sigma) };
            const double eps{ optimizer.eps() };

            for (std::size_t i = 0; i < numStrikes; ++i)
            {
                // Vector of calendar constraints context
                const double k{ kSlice[i] };
                const double wPrev{ detail::calculateWk(a0, b0, rho0, m0, sigma0, k) };
                calCtx[i] = { k, wPrev, eps, atmWK };
            }

            addCalendarConstraint(optimizer, calCtx);
        }

        // ---------- Convexity constraints ----------

        Vector<ConvexityContexts> convexityCtxs;
        convexityCtxs.reserve(numStrikes);
        for (double k : kSlice)
        {
            convexityCtxs.push_back({ k, atmWK });
            addConvexityConstraint(optimizer, convexityCtxs.back());
        }

        // ---------- Set objective function ----------

        ObjectiveContexts obj
        {
            kSlice.data(),
            wKSlice.data(),
            numStrikes,
            optimizer.userValue()
        };

        setMinObjective(optimizer, obj);

        // ---------- Run optimization ----------

        Vector<double> params{ optimizer.optimize() };

        return Params{
            T,
            aParam(atmWK,
                   Real(params[0]),
                   Real(params[1]),
                   Real(params[2]),
                   Real(params[3])),
            Real(params[0]),
            Real(params[1]),
            Real(params[2]),
            Real(params[3])
        };
    }

    template <::nlopt::algorithm Algo>
    void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
    {
        optimizer.addInequalityConstraint(
            +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

                // Extract data
                const double eps{ opt.eps() };
                const double atmWK{ opt.userValue() };
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                // Precomputes
                const double R{ std::sqrt(m * m + sigma * sigma) };
                const double S{ std::sqrt(1.0 - rho * rho) };
              
                // Calculate
                const double a{ aParam(atmWK, b, rho, m, sigma) };
                const double wMin{ std::fma(b, sigma * S, a) };

                if (grad)
                {
                    grad[0] = -(rho * m - R + sigma * S);       
                    grad[1] = -(b * (m - sigma * rho / S));      
                    grad[2] = -(b * (rho - m / R));            
                    grad[3] = -(b * (S - sigma / R));            
                }

                return eps - wMin;
            },
            &optimizer
        );
    }

    template <::nlopt::algorithm Algo>
    void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
    {
        // ---------- Right wing: b*(1 + rho) >= eps ----------
        
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void* data) -> double
            {  
                // Reinterpret opaque pointer
                const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

                // Extract data
                const double eps{ opt.eps() };
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = -(1.0 + rho);  
                    grad[1] = -b;             
                    grad[2] = 0.0;
                    grad[3] = 0.0;
                }

                return eps - b * (1.0 + rho);
            },
            &optimizer
        );

        // ---------- Left wing: b*(1 - rho) >= eps ----------

        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

                // Extract data
                const double eps{ opt.eps() };
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = -(1.0 - rho);   
                    grad[1] = b;              
                    grad[2] = 0.0;
                    grad[3] = 0.0;
                }

                return eps - b * (1.0 - rho);
            },
            &optimizer
        );
    }

    template <::nlopt::algorithm Algo>
    void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
    {
        // ---------- Right wing: b*(1 + rho) <= 2 ----------

        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                // Extract data
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = (1.0 + rho);   
                    grad[1] = b;            
                    grad[2] = 0.0;           
                    grad[3] = 0.0;           
                }

                return b * (1.0 + rho) - 2.0;
            },
            nullptr
        );

        // ---------- Left wing: b*(1 - rho) <= 2 ----------

        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                // Extract data
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = (1.0 - rho);   
                    grad[1] = -b;            
                    grad[2] = 0.0;          
                    grad[3] = 0.0;           
                }

                return b * (1.0 - rho) - 2.0;
            },
            nullptr
        );
    }
  

    template <::nlopt::algorithm Algo>
    void addCalendarConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer,
        Vector<CalendarContexts>& ctx)
    {
        for (auto& ctxK : ctx)
        {
            optimizer.addInequalityConstraint(
                +[](unsigned, const double* x, double* grad, void* data) -> double
                {
                    // Reinterpret opaque pointer
                    const CalendarContexts& ctxK{ *static_cast<const CalendarContexts*>(data) };

                    // Extract data
                    const double atmWK{ ctxK.atmWK };
                    const double k{ ctxK.k };
                    const double prevWk{ ctxK.prevWk };
                    const double eps{ ctxK.eps };
                    const double b{ x[0] };
                    const double rho{ x[1] };
                    const double m{ x[2] };
                    const double sigma{ x[3] };

                    // Precomputes
                    const double xi{ k - m };
                    const double Rk{ std::hypot(xi, sigma) };
                    const double invRk{ 1.0 / Rk };
                    const double R0{ std::sqrt(m * m + sigma * sigma) };
                    const double invR0{ 1.0 / R0 };

                    // Calculate
                    const double a{ aParam(atmWK, b, rho, m, sigma) };
                    const double wK{ calculateWk(a, b, rho, m, sigma, k) };

                    if (grad)
                    {
                        grad[0] = -(rho * xi + Rk - (-rho * m + R0));        
                        grad[1] = -b * (xi + m);                             
                        grad[2] = b * (m * invR0 + xi * invRk);              
                        grad[3] = -b * (sigma * invRk - sigma * invR0);      
                    }

                    return prevWk + eps - wK;
                },
                &ctxK
            );
        }
    }

    template <::nlopt::algorithm Algo>
    void addConvexityConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer, 
        ConvexityContexts& ctx)
    {
        optimizer.addInequalityConstraint(
            +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const auto& ctx = *static_cast<const ConvexityContexts*>(data);

                // Extract data
                const double k{ ctx.k };
                const double atmWK{ ctx.atmWK };
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                // Calculate
                const double a{ aParam(atmWK, b, rho, m, sigma) };
                const GkCache p{ a, b, rho, m, sigma, k };
                const double g{ gk(p) };

                if (grad)
                {
                    const auto dg{ gkGrad(b, rho, m, sigma, k, p) };

                    grad[0] = -dg[0]; 
                    grad[1] = -dg[1]; 
                    grad[2] = -dg[2]; 
                    grad[3] = -dg[3]; 
                }

                return -g;
            },
            &ctx
        );
    }

    template <::nlopt::algorithm Algo>
    void setMinObjective(opt::nlopt::Optimizer<4, Algo>& optimizer,  
        const ObjectiveContexts& ctx)
    {
        optimizer.setMinObjective(
            +[](unsigned n, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const ObjectiveContexts& ctx = *static_cast<const ObjectiveContexts*>(data);

                // Extract data
                const double atmWK{ ctx.atmWK };
                std::size_t N{ ctx.n };
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                // Precomputes
                const double R0{ std::sqrt(m * m + sigma * sigma) };
                const double invR0{ 1.0 / R0 };

                // Calculate
                const double a{ aParam(atmWK, b, rho, m, sigma) };

                // Initialize
                double SSE{ 0.0 };
                std::array<double, 4> g{}; 

                for (std::size_t i = 0; i < N; ++i)
                {
                    // Extract data
                    const double k{ ctx.k[i] };
                    const double wM{ ctx.wM[i] };

                    // Precomputes
                    const double xi{ k - m };
                    const double R{ std::hypot(xi, sigma) };

                    // Calculate
                    const double wK{ calculateWk(a, b, rho, m, sigma, k) };
                    const double r{ wK - wM };
                    SSE += r * r;

                    if (grad)
                    {
                        // Precomputes
                        const double invR{ 1.0 / R };

                        g[0] += 2.0 * r * ((rho * xi + R) - (-rho * m + R0));
                        g[1] += 2.0 * r * (b * (xi + m));
                        g[2] += 2.0 * r * (-b * (m * invR0 + xi * invR));
                        g[3] += 2.0 * r * (b * (sigma * invR - sigma * invR0));
                    }
                }

                if (grad)
                {
                    const unsigned mCopy = (std::min)(n, static_cast<unsigned>(g.size()));
                    std::copy_n(g.data(), mCopy, grad);
                }

                return SSE;
            },
            const_cast<ObjectiveContexts*>(&ctx)
        );
    }

    template <std::floating_point T>
    T calculateWk(T a,
        T b,
        T rho,
        T m,
        T sigma,
        T k) noexcept
    {
        const T x{ k - m };
        return std::fma(b, (rho * x + std::hypot(x, sigma)), a);
    }
} 