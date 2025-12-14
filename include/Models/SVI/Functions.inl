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


#include "Utils/IO/Log.hpp"

#include <algorithm>   
#include <array>       
#include <cmath>      
#include <cstddef>     
#include <format>     
#include <limits>    
#include <utility>     


namespace uv::models::svi
{
    template <::nlopt::algorithm Algo>
    std::tuple<std::vector<Params>, core::VolSurface> calibrate(const core::VolSurface& mktVolSurf,
        const opt::nlopt::Optimizer<4, Algo>& prototype,
        bool isValidateResults)
    {
        // Copy market volatility surface
        core::VolSurface sviVolSurf{ mktVolSurf };

        // Initialize vectors
        std::vector<Params> sviSlices;
        sviSlices.reserve(sviVolSurf.numTenors());

        bool isFirstSlice{ true };

        // Calibrate each slice
        for (auto& slice : sviVolSurf.slices())
        {
            // Extract and convert data
            const Vector<double> wKSlice{slice.wT().cbegin(), slice.wT().cend() };
            const Vector<double> kSlice{slice.logKF().cbegin(), slice.logKF().cend()};

            // Initialize optimizer instance
            opt::nlopt::Optimizer optimizer{ prototype.fresh() };
            optimizer.setUserValue(static_cast<double>(slice.atmWT()));
            

            // Set initial guess and bounds
            optimizer.setGuessBounds
            (
                detail::initGuess(slice),
                detail::lowerBounds(slice),
                detail::upperBounds(slice)
            );

            // Enforce positive minimum total variance constraint
            detail::addWMinConstraint(optimizer);

            // Enforce Roger Lee wing slope constraints
            detail::addMinSlopeConstraint(optimizer);
            detail::addMaxSlopeConstraint(optimizer);

            std::vector<detail::CalendarContexts> calCtx;

            if (!isFirstSlice)
            {
                calCtx.resize(kSlice.size());
                const Params& prevParams{ sviSlices.back() };

                const double a0 = double(prevParams.a);
                const double b0 = double(prevParams.b);
                const double rho0 = double(prevParams.rho);
                const double m0 = double(prevParams.m);
                const double sigma0 = double(prevParams.sigma);

                for (std::size_t i = 0; i < kSlice.size(); ++i)
                {
                    const double k = kSlice[i];
                    const double xi = k - m0;
                    const double R = std::hypot(xi, sigma0);
                    const double wPrev = a0 + b0 * (rho0 * xi + R);

                    calCtx[i] = { k, wPrev, optimizer.eps(), optimizer.userValue() };
                }

                detail::addCalendarConstraint(optimizer, calCtx);
            }

            // Enforce convexity constraints: g(k) ≥ 0
            std::vector<detail::ConvexityContexts> convexityCtxs;
            convexityCtxs.reserve(kSlice.size());
            for (double k : kSlice)
            {
                convexityCtxs.push_back(
                    detail::ConvexityContexts{
                        k,
                        optimizer.userValue()   // ATM total variance
                    }
                );

                detail::addConvexityConstraint(
                    optimizer,
                    convexityCtxs.back()
                );
            }

            // Objective function contexts
            detail::ObjectiveContexts obj{ kSlice.data(), wKSlice.data(), kSlice.size(), optimizer.userValue() };

            // Define objective function with analytical gradient
            detail::setMinObjective(optimizer, obj);

            // Solve the optimization problem
            Vector<double> params{ optimizer.optimize() };

            // Extract calibration results
            Real T{ slice.T() };
            Real b{ Real(params[0]) };
            Real rho{ Real(params[1]) };
            Real m{ Real(params[2]) };
            Real sigma{ Real(params[3]) };
            const Real a
            {
                slice.atmWT() - b * (-rho * m + std::sqrt(m * m + sigma * sigma))
            };

            Params sviSlice{ T, a, b, rho, m, sigma };

            // Evaluate calibration
            //if (isValidateResults) detail::evalCal(sviSlice, optimizer, kSlice, wkSlice);

            // Save calibration parameter results
            sviSlices.emplace_back(std::move(sviSlice));

            // Update calculated variances use them on the next slice calibration
            Vector<double> wkSlice = detail::makewKSlice<double>(kSlice, double(a), double(b), double(rho), double(m), double(sigma));

            // Make a local copy first
            Vector<Real> wtCopy{ wkSlice.begin(), wkSlice.end() };

            // Set the slice of the calibrated volatility surface
            slice.setWT(wtCopy);

            isFirstSlice = false;
        }
        return { std::move(sviSlices), std::move(sviVolSurf) };
    }
} // namespace uv::models::svi

namespace uv::models::svi::detail
{
    struct ObjectiveContexts
    {
        const double* k;      // Pointer to log-forward moneyness
        const double* wM;     // Pointer to marekt total variance
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
        double atmWK;   // Atm varaince
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
    void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
    {
        optimizer.addInequalityConstraint(
            +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

                // Extract variables
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

                // Extract variables
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

                // Extract variables
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
                // Extract variables
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
                // Extract variables
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
        std::vector<CalendarContexts>& ctx)
    {
        for (auto& ctxK : ctx)
        {
            optimizer.addInequalityConstraint(
                +[](unsigned, const double* x, double* grad, void* data) -> double
                {
                    // Reinterpret opaque pointer
                    const CalendarContexts& ctxK{ *static_cast<const CalendarContexts*>(data) };

                    // Extract variables
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

                // Extract variables
                const double k{ ctx.k };
                const double atmWK{ ctx.atmWK };
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                // Precomputes
                const double R0{ std::sqrt(m * m + sigma * sigma) };

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

                // Extract variables
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
                    // Extract variables
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


    template <::nlopt::algorithm Algo>
    void evalCal(const Params& sviSlice,
        const opt::nlopt::Optimizer<4, Algo>& optimizer,
        const Vector<double>& kSlice,
        const Vector<double>& wKPrevSlice) noexcept
    {
        // Extract parameters
        const double a{ double(sviSlice.a) };
        const double b{ double(sviSlice.b) };
        const double rho{ double(sviSlice.rho) };
        const double m{ double(sviSlice.m) };
        const double sigma{ double(sviSlice.sigma) };

        // Extract config attributes
        const double eps{ optimizer.eps() };
        const double tol{ optimizer.tol() };

        // ---- wMin constraint violation ----
        const double S = std::sqrt(std::max(0.0, 1.0 - rho * rho));
        const double wMin = std::fma(b, sigma * S, a);

        UV_WARN(eps - wMin > tol,
            std::format("No-arbitrage: wMin violated: wMin = {:.6f} < eps = {:.6f}",
                wMin, eps));

        // ---- Lee wing-slope constraint violations ----
        const double leeRight = b * (1.0 + rho) - 2.0;
        const double leeLeft = b * (1.0 - rho) - 2.0;

        UV_WARN(leeRight > tol,
            std::format("No-arbitrage: right wing slope > 2: b*(1+rho) = {:.6f} "
                "(b = {:.6f}, rho = {:.6f})",
                b * (1.0 + rho), b, rho));

        UV_WARN(leeLeft > tol,
            std::format("No-arbitrage: left wing slope > 2: b*(1-rho) = {:.6f} "
                "(b = {:.6f}, rho = {:.6f})",
                b * (1.0 - rho), b, rho));

        // ---- Convexity g(k) violation ----
        double gMin = std::numeric_limits<double>::infinity();
        double kAtMin = 0.0;
        std::size_t nViol = 0;

        for (double k : kSlice)
        {
            const GkCache pre{ a, b, rho, m, sigma, k };
            const double g = gk(pre);
            if (g < gMin) { gMin = g; kAtMin = k; }
            if (g < -tol) ++nViol;
        }

        UV_WARN(gMin < -tol,
            std::format("No-arbitrage: convexity violated: min g(k) = {:.6e} at k = {:.4f} "
                "(violations {} / {}, tol = {:.2e})",
                gMin, kAtMin, nViol, kSlice.size(), tol));

        //// ---- Calendar no-arb: w_curr(k) >= w_prev(k) + eps ----
        //std::size_t nCalViol = 0;
        //double minMargin = std::numeric_limits<double>::infinity();
        //double kAtWorst = 0.0;

        //for (std::size_t i = 0; i < kSlice.size(); ++i)
        //{
        //    const double k = kSlice[i];
        //    const double wPrev = wKPrevSlice[i] + eps;
        //    const double wCurr = wk(a, b, rho, m, sigma, k);
        //    const double margin = wCurr - wPrev; // should be >= 0
        //    if (margin < minMargin) { minMargin = margin; kAtWorst = k; }
        //    if (margin < -tol) ++nCalViol;
        //}

        //UV_WARN(minMargin < -tol,
        //    std::format("No-arbitrage: calendar violated: "
        //        "min[w_curr - (w_prev+eps)] = {:.6e} at k = {:.4f} "
        //        "(violations {} / {}, tol = {:.2e}, eps = {:.2e})",
        //        minMargin, kAtWorst, nCalViol, kSlice.size(),
        //        tol, eps));
    }

    template<std::floating_point T>
    Vector<T> makewKSlice(const Vector<T>& kSlice,
        T a, T b, T rho, T m, T sigma) noexcept
    {
        Vector<T> wKSlice;
        wKSlice.reserve(kSlice.size());

        std::transform(
            kSlice.begin(), kSlice.end(),
            std::back_inserter(wKSlice),
            [a, b, rho, m, sigma](T k) noexcept
            {
                return calculateWk(a, b, rho, m, sigma, k);
            });

        return wKSlice;
    }
} 
