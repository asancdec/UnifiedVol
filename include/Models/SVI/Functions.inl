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

#include <iostream>

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
        Vector<double> wkSlice(sviVolSurf.numStrikes(), 0.0);

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
            detail::addMinSlopeConstraint(optimizer, 1e-4);
            detail::addMaxSlopeConstraint(optimizer);

            if (!isFirstSlice)
            {
                // Vector of constraints context
                std::vector<detail::CalendarCtx> contexts;
                contexts.resize(kSlice.size());

                Params& prevParams{ sviSlices.back()};

                // Evaluate
                const double a{ double(prevParams.a) };
                const double b{ double(prevParams.b) };
                const double rho{ double(prevParams.rho) };
                const double m{ double(prevParams.m) };
                const double sigma{ double(prevParams.sigma) };


                for (std::size_t i = 0; i < kSlice.size(); ++i)
                {
                    // Extract preivous params
                    const double k{ kSlice[i] };
                    const double xi{ k - m };
                    const double R{ std::hypot(xi, sigma) };

                    const double wKPrev
                    {
                        a + b * (rho * xi + R)
                    };
                    contexts[i] = detail::CalendarCtx
                    {
                        k,
                        wKPrev,
                        optimizer.eps(),
                        optimizer.userValue()
                    };
                }

                // Enforce calendar spread arbitrage constraints: Wk_current ≥ Wk_previous
                detail::addCalendarConstraint(optimizer, contexts);
            }

            // Enforce convexity constraints: g(k) ≥ 0
            std::vector<detail::ConvexityCtx> convexityCtxs;
            convexityCtxs.reserve(kSlice.size());
            for (double k : kSlice)
            {
                convexityCtxs.push_back(
                    detail::ConvexityCtx{
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
            detail::ObjCtx obj{ kSlice.data(), wKSlice.data(), kSlice.size(), optimizer.userValue() };

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

            std::cout << "alpha" << a << "\n";

            Params sviSlice{ T, a, b, rho, m, sigma };

            // Evaluate calibration
            if (isValidateResults) detail::evalCal(sviSlice, optimizer, kSlice, wkSlice);

            // Save calibration parameter results
            sviSlices.emplace_back(std::move(sviSlice));

            // Update calculated variances use them on the next slice calibration
            wkSlice = detail::makewKSlice<double>(kSlice, double(a), double(b), double(rho), double(m), double(sigma));

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
    struct ObjCtx
    {
        const double* k;     // Pointer to log-forward moneyness
        const double* wK;    // Pointer to total variance
        std::size_t   n;     // Number of points
        const double atmWT;
    };

    struct CalendarCtx
    {
        double k;       // Log-forward moneyness
        double prevWk;  // Total variance of the previous slice
        double eps;     // Epsilon value
        double atmWT;   // Atm vol
    };

    struct ConvexityCtx
    {
        double k;
        double atmWT;   // Atm vol
    };

    struct GKPrecomp
    {
        double x;             // x := k-m
        double R;             // R:= sqrt(x^2 + sigma^2)
        double invR;          // invR := 1 / R
        double wk;            // w(k) = a + b*(rho*x + R)
        double wkD1;          // w'(k) = b * (rho + x/R)
        double wkD1Squared;   // w'(k)^2
        double invRCubed;     // 1/(R^3)
        double sigmaSquared;  // sigma^2
        double wkD2;          // w''(k) = b * sigma^2 / R^3
        double A;             // A := 1 - k * w'/(2 * w)                        
        double B;             // B := 1/w(k) + 1/4

        explicit GKPrecomp(double a, double b, double rho, double m, double sigma, double k) noexcept
        {
            this->x = k - m;                                      // x := k-m
            this->R = std::hypot(x, sigma);                      // R:= sqrt(x^2 + sigma^2)
            this->invR = 1.0 / R;                                // invR := 1 / R
            this->wk = std::fma(b, (rho * x + R), a);            // w(k) = a + b*(rho*x + R)
            this->wkD1 = b * (rho + x * invR);                   // w'(k) = b * (rho + x/R)
            this->wkD1Squared = wkD1 * wkD1;                     // w'(k)^2
            this->invRCubed = invR * invR * invR;                // 1/(R^3)
            this->sigmaSquared = sigma * sigma;                  // sigma^2
            this->wkD2 = b * sigmaSquared * invRCubed;           // w''(k) = b * sigma^2 / R^3
            this->A = 1.0 - 0.5 * k * wkD1 / wk;                 // A := 1 - k * w'/(2 * w)                        
            this->B = (1.0 / wk) + 0.25;                         // B := 1/w(k) + 1/4
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

                const double eps{ opt.eps() };
                const double wATM{ opt.userValue() };

                // Extract variables
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                // Precompute
                const double R{ std::sqrt(m * m + sigma * sigma) };
                const double S{ std::sqrt(1.0 - rho * rho) };

                const double a{wATM - b * (-rho * m + R)};

                // w_min = a + b * sigma * sqrt(1 - rho^2)
                const double wMin{ std::fma(b, sigma * S, a) };

                if (grad)
                {
                    grad[0] = -(rho * m - R + sigma * S);        // ∂/∂b
                    grad[1] = -(b * (m - sigma * rho / S));      // ∂/∂rho
                    grad[2] = -(b * (rho - m / R));              // ∂/∂m
                    grad[3] = -(b * (S - sigma / R));            // ∂/∂sigma
                }

                // Enforce wMin ≥ eps  c(x) = eps - wMin ≤ 0
                return eps - wMin;
            },
            &optimizer
        );
    }

    template <::nlopt::algorithm Algo>
    void addMinSlopeConstraint(
        opt::nlopt::Optimizer<4, Algo>& optimizer,
        double epsSlope
    )
    {
        // Right wing: b*(1 + rho) >= epsSlope
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void* data) -> double
            {
                const double eps = *static_cast<const double*>(data);

                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = -(1.0 + rho);   // ∂/∂b
                    grad[1] = -b;             // ∂/∂rho
                    grad[2] = 0.0;
                    grad[3] = 0.0;
                }

                // c(x) = epsSlope - b*(1+rho) <= 0
                return eps - b * (1.0 + rho);
            },
            &epsSlope
        );

        // Left wing: b*(1 - rho) >= epsSlope
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void* data) -> double
            {
                const double eps = *static_cast<const double*>(data);

                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = -(1.0 - rho);   // ∂/∂b
                    grad[1] = b;              // ∂/∂rho
                    grad[2] = 0.0;
                    grad[3] = 0.0;
                }

                // c(x) = epsSlope - b*(1-rho) <= 0
                return eps - b * (1.0 - rho);
            },
            &epsSlope
        );
    }

    template <::nlopt::algorithm Algo>
    void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
    {
        // Right wing: b*(1 + rho) <= 2
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                // x = [ b, rho, m, sigma ]
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = (1.0 + rho);   // ∂/∂b
                    grad[1] = b;             // ∂/∂rho
                    grad[2] = 0.0;           // ∂/∂m
                    grad[3] = 0.0;           // ∂/∂sigma
                }

                // c(x) = b * (1 + rho) - 2 <= 0
                return b * (1.0 + rho) - 2.0;
            },
            nullptr
        );

        // Left wing: b*(1 - rho) <= 2
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                // x = [ b, rho, m, sigma ]
                const double b{ x[0] };
                const double rho{ x[1] };

                if (grad)
                {
                    grad[0] = (1.0 - rho);   // ∂/∂b
                    grad[1] = -b;            // ∂/∂rho
                    grad[2] = 0.0;           // ∂/∂m
                    grad[3] = 0.0;           // ∂/∂sigma
                }

                // c(x) = b * (1 - rho) - 2 <= 0
                return b * (1.0 - rho) - 2.0;
            },
            nullptr
        );
    }
  

    template <::nlopt::algorithm Algo>
    void addCalendarConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer,
        std::vector<CalendarCtx>& contexts)
    {
        // For each k, add one no calendar arbitrage constraint
        for (auto& ctx : contexts)
        {
            optimizer.addInequalityConstraint(
                +[](unsigned, const double* x, double* grad, void* data) -> double
                {
                    // Reinterpret opaque pointer
                    const CalendarCtx& c{ *static_cast<const CalendarCtx*>(data) };

                    // Extract variables
                    const double b{ x[0] };
                    const double rho{ x[1] };
                    const double m{ x[2] };
                    const double sigma{ x[3] };

                    // Precompute variables
                    const double xi{ c.k - m };
                    const double Rk{ std::hypot(xi, sigma) };
                    const double invRk{ 1.0 / Rk };

                    const double R0{ std::sqrt(m * m + sigma * sigma) };
                    const double invR0{ 1.0 / R0 };

                    const double a{c.atmWT - b * (-rho * m + R0)};


                    // Calculate wK
                    const double wK{ a + b * (rho * xi + Rk) };

                    if (grad)
                    {
                        grad[0] = -(rho * xi + Rk - (-rho * m + R0));        // ∂/∂b
                        grad[1] = -b * (xi + m);                             // ∂/∂rho
                        grad[2] = b * (m * invR0 + xi * invRk);              // ∂/∂m
                        grad[3] = -b * (sigma * invRk - sigma * invR0);      // ∂/∂sigma
                    }

                    // NLopt expects c(x) ≤ 0.
                    // Here:  c(x) = prevWk + eps − w(k)
                    // Enforces w(k) ≥ prevWk + eps  (no calendar arbitrage, with small tolerance)
                    return c.prevWk + c.eps - wK;
                },
                &ctx
            );
        }
    }

    template <::nlopt::algorithm Algo>
    void addConvexityConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer, ConvexityCtx& ctx)
    {
        optimizer.addInequalityConstraint(
            +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
            {
                //--------------------------------------------------------------------------
                // Context
                //--------------------------------------------------------------------------
                const auto& c =
                    *static_cast<const ConvexityCtx*>(data);

                const double k{ c.k };
                const double wATM{ c.atmWT };

                //--------------------------------------------------------------------------
                // Extract variables
                //   x = [ b, rho, m, sigma ]
                //--------------------------------------------------------------------------
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                //--------------------------------------------------------------------------
                // Recompute a from ATM
                //--------------------------------------------------------------------------
                const double R0{ std::sqrt(m * m + sigma * sigma) };
                const double a
                {
                    wATM - b * (-rho * m + R0)
                };

                //--------------------------------------------------------------------------
                // Convexity g(k)
                //--------------------------------------------------------------------------
                const GKPrecomp p{ a, b, rho, m, sigma, k };
                const double g{ gk(p) };

                if (grad)
                {
                    const auto dg{ gkGrad(a, b, rho, m, sigma, k, p) };

                    // c(x) = -g(k)
                    grad[0] = -dg[0]; // ∂/∂b
                    grad[1] = -dg[1]; // ∂/∂rho
                    grad[2] = -dg[2]; // ∂/∂m
                    grad[3] = -dg[3]; // ∂/∂sigma
                }

                // Enforce g(k) ≥ 0
                return -g;
            },
            &ctx
        );
    }

    template <::nlopt::algorithm Algo>
    void setMinObjective(opt::nlopt::Optimizer<4, Algo>& optimizer,  const ObjCtx& obj) 
    {
        optimizer.setMinObjective(
            +[](unsigned n, const double* x, double* grad, void* data) -> double
            {
                //--------------------------------------------------------------------------
                // Context
                //--------------------------------------------------------------------------
                const ObjCtx& obj =
                    *static_cast<const ObjCtx*>(data);

                const double wATM{ obj.atmWT };

                //--------------------------------------------------------------------------
                // Extract parameters
                //   x = [ b, rho, m, sigma ]
                //--------------------------------------------------------------------------
                const double b{ x[0] };
                const double rho{ x[1] };
                const double m{ x[2] };
                const double sigma{ x[3] };

                //--------------------------------------------------------------------------
                // Reconstruct a from ATM
                //--------------------------------------------------------------------------
                const double R0{ std::sqrt(m * m + sigma * sigma) };
                const double a
                {
                    wATM - b * (-rho * m + R0)
                };

                //--------------------------------------------------------------------------
                // Objective and gradient accumulator
                //--------------------------------------------------------------------------
                double SSE{ 0.0 };
                std::array<double, 4> g{}; // [b, rho, m, sigma]

                //--------------------------------------------------------------------------
                // Loop over strikes
                //--------------------------------------------------------------------------
                for (std::size_t i = 0; i < obj.n; ++i)
                {
                    const double k{ obj.k[i] };
                    const double wM{ obj.wK[i] };

                    const double xi{ k - m };
                    const double R{ std::hypot(xi, sigma) };

                    const double wKModel
                    {
                        a + b * (rho * xi + R)
                    };

                    const double r{ wKModel - wM };
                    SSE += r * r;

                    if (grad)
                    {
                        const double invR{ 1.0 / R };
                        const double invR0{ 1.0 / R0 };

                        // Derivatives of w(k)
                        const double dw_db =
                            (rho * xi + R) - (-rho * m + R0);

                        const double dw_drho =
                            b * (xi + m);

                        const double dw_dm =
                            -b * (m * invR0 + xi * invR);

                        const double dw_dsigma =
                            b * (sigma * invR - sigma * invR0);

                        g[0] += 2.0 * r * dw_db;
                        g[1] += 2.0 * r * dw_drho;
                        g[2] += 2.0 * r * dw_dm;
                        g[3] += 2.0 * r * dw_dsigma;
                    }
                }

                //--------------------------------------------------------------------------
                // Write gradient
                //--------------------------------------------------------------------------
                if (grad)
                {
                    const unsigned mCopy =
                        (std::min)(n, static_cast<unsigned>(g.size()));

                    std::copy_n(g.data(), mCopy, grad);
                }

                return SSE;
            },
            const_cast<ObjCtx*>(&obj)
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
            const GKPrecomp pre{ a, b, rho, m, sigma, k };
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
                return wk(a, b, rho, m, sigma, k);
            });

        return wKSlice;
    }
} 
