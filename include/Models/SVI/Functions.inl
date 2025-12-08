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
        const opt::nlopt::Optimizer<5, Algo>& prototype,
        bool isValidateResults)
    {
        // Copy market volatility surface
        core::VolSurface sviVolSurf{ mktVolSurf };

        // Initialize vectors
        std::vector<Params> sviSlices;
        sviSlices.reserve(sviVolSurf.numTenors());
        Vector<double> wkSlice(sviVolSurf.numStrikes(), 0.0);

        // Calibrate each slice
        for (auto& slice : sviVolSurf.slices())
        {
            // Extract and convert data
            const Vector<double> wKSlice{slice.wT().cbegin(), slice.wT().cend() };
            const Vector<double> kSlice{slice.logKF().cbegin(), slice.logKF().cend()};

            // Initialize optimizer instance
            opt::nlopt::Optimizer optimizer{ prototype.fresh() };

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
            detail::addMaxSlopeConstraint(optimizer);

            // Vector of constraints context
            std::vector<detail::ConstraintCtx> contexts;
            contexts.resize(kSlice.size());
            for (std::size_t i = 0; i < kSlice.size(); ++i)
            {
                contexts[i] = detail::ConstraintCtx
                {
                    kSlice[i],
                    wkSlice[i],
                    optimizer.eps()
                };
            }

            // Enforce calendar spread arbitrage constraints: Wk_current ≥ Wk_previous
            detail::addCalendarConstraint(optimizer, contexts);

            // Enforce convexity constraints: g(k) ≥ 0
            Vector<double> kStorage{ kSlice };
            detail::addConvexityConstraint(optimizer, kStorage);

            // Objective function contexts
            detail::ObjCtx obj{ kSlice.data(), wKSlice.data(), kSlice.size() };

            // Define objective function with analytical gradient
            detail::setMinObjective(optimizer, obj);

            // Solve the optimization problem
            Vector<double> params{ optimizer.optimize() };

            // Extract calibration results
            Real T{ slice.T() };
            Real a{ Real(params[0]) };
            Real b{ Real(params[1]) };
            Real rho{ Real(params[2]) };
            Real m{ Real(params[3]) };
            Real sigma{ Real(params[4]) };

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
    };

    struct ConstraintCtx
    {
        double k;       // Log-forward moneyness
        double prevWk;  // Total variance of the previous slice
        double eps;     // Epsilon value
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
    void addWMinConstraint(opt::nlopt::Optimizer<5, Algo>& optimizer) noexcept
    {
        const double* epsPtr{ &optimizer.eps() };

        optimizer.addInequalityConstraint(
            +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const double eps{ *static_cast<const double*>(data) };

                // Extract variables
                const double a{ x[0] };
                const double b{ x[1] };
                const double rho{ x[2] };
                const double sigma{ x[4] };

                // Precompute
                const double S{ std::sqrt(1.0 - rho * rho) };

                // w_min = a + b * sigma * sqrt(1 - rho^2)
                const double wMin{ std::fma(b, sigma * S, a) };

                if (grad)
                {
                    // c(x) = eps - wMin  =>  ∇c = -∇wMin
                    grad[0] = -1.0;                         // -∂wMin/∂a
                    grad[1] = -sigma * S;                   // -∂wMin/∂b
                    grad[2] = -b * sigma * (rho / S);       // -∂wMin/∂rho (since ∂S/∂rho = -rho/S)
                    grad[3] = 0.0;                          // -∂wMin/∂m
                    grad[4] = -b * S;                       // -∂wMin/∂sigma
                }

                // Enforce wMin ≥ eps  c(x) = eps - wMin ≤ 0
                return eps - wMin;
            },
            const_cast<double*>(epsPtr) // nlopt API takes void*
        );
    }

    template <::nlopt::algorithm Algo>
    void addMaxSlopeConstraint(opt::nlopt::Optimizer<5, Algo>& optimizer) noexcept
    {
        // Right wing: b*(1 + rho) <= 2
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                // x = [a, b, rho, m, sigma]
                const double b{ x[1] };
                const double rho{ x[2] };

                if (grad)
                {
                    grad[0] = 0.0;            // d/da
                    grad[1] = (1.0 + rho);    // d/db
                    grad[2] = b;              // d/drho
                    grad[3] = 0.0;            // d/dm
                    grad[4] = 0.0;            // d/dsigma
                }

                // c(x) = b*(1 + rho) - 2 <= 0
                return b * (1.0 + rho) - 2.0;
            },
            nullptr
        );

        // Left wing: b*(1 - rho) <= 2
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void*) -> double
            {
                const double b{ x[1] };
                const double rho{ x[2] };

                if (grad)
                {
                    grad[0] = 0.0;            // d/da
                    grad[1] = (1.0 - rho);    // d/db
                    grad[2] = -b;             // d/drho
                    grad[3] = 0.0;            // d/dm
                    grad[4] = 0.0;            // d/dsigma
                }

                // c(x) = b*(1 - rho) - 2 <= 0
                return b * (1.0 - rho) - 2.0;
            },
            nullptr
        );
    }

    template <::nlopt::algorithm Algo>
    void addCalendarConstraint(opt::nlopt::Optimizer<5, Algo>& optimizer,
        std::vector<ConstraintCtx>& contexts) noexcept
    {
        // For each k, add one no calendar arbitrage constraint
        for (auto& ctx : contexts)
        {
            optimizer.addInequalityConstraint(
                +[](unsigned, const double* x, double* grad, void* data) -> double
                {
                    // Reinterpret opaque pointer
                    const ConstraintCtx& c{ *static_cast<const ConstraintCtx*>(data) };

                    // Extract variables
                    const double a{ x[0] };
                    const double b{ x[1] };
                    const double rho{ x[2] };
                    const double m{ x[3] };
                    const double sigma{ x[4] };

                    // Precompute variables
                    const double xi{ c.k - m };
                    const double R{ std::hypot(xi, sigma) };
                    const double invR{ 1.0 / R };

                    // Calculate wK
                    const double wK{ a + b * (rho * xi + R) };

                    if (grad)
                    {
                        grad[0] = -1.0;                     // -∂w/∂a
                        grad[1] = -(rho * xi + R);          // -∂w/∂b
                        grad[2] = -b * xi;                  // -∂w/∂ρ
                        grad[3] = b * (rho + xi * invR);    // -∂w/∂m
                        grad[4] = -b * (sigma * invR);      // -∂w/∂σ
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
    void addConvexityConstraint(opt::nlopt::Optimizer<5, Algo>& optimizer, Vector<double>& kStorage) noexcept
    {
        // For each k, add one convexity constraint g(k) ≥ 0
        for (std::size_t i = 0; i < kStorage.size(); ++i)
        {
            optimizer.addInequalityConstraint(
                +[](unsigned n, const double* x, double* grad, void* data) -> double
                {
                    // Reinterpret opaque pointer
                    const double k{ *static_cast<const double*>(data) };

                    // Extract variables
                    const double a{ x[0] };
                    const double b{ x[1] };
                    const double rho{ x[2] };
                    const double m{ x[3] };
                    const double sigma{ x[4] };

                    // Build g(K) precompute instance for efficiency
                    const GKPrecomp p{ a, b, rho, m, sigma, k };
                    const double g{ gk(p) };

                    if (grad)
                    {
                        // ∇g = [dg/da, dg/db, dg/dρ, dg/dm, dg/dσ]
                        const auto dg{ gkGrad(a, b, rho, m, sigma, k, p) };
                        for (unsigned j = 0; j < 5 && j < n; ++j)
                            grad[j] = -dg[j]; // c(x) = -g(k)
                    }

                    // NLopt enforces c(x) ≤ tol. 
                    // Here c(x) = -g(k), so -g(k) ≤ tol  
                    // Enforces g(k) ≥ -tol  (≈ g(k) ≥ 0 with small slack)
                    return -g;
                },
                &kStorage[i]  // pointer to stable local copy
            );
        }
    }

    template <::nlopt::algorithm Algo>
    void setMinObjective(opt::nlopt::Optimizer<5, Algo>& optimizer, const ObjCtx& obj) noexcept
    {
        optimizer.setMinObjective(
            +[](unsigned n, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const ObjCtx& obj{ *static_cast<const ObjCtx*>(data) };

                // Extract variables
                const double a{ x[0] };
                const double b{ x[1] };
                const double rho{ x[2] };
                const double m{ x[3] };
                const double sigma{ x[4] };

                // Initialize variable to store total SSE
                double SSE{ 0.0 };

                // Accumulate ∇f here
                std::array<double, 5> g{};

                for (std::size_t i = 0; i < obj.n; ++i)
                {
                    // Extract variables
                    const double k{ obj.k[i] };
                    const double wKMarket{ obj.wK[i] };

                    // Precompute variables
                    const double xi{ k - m };
                    const double R{ std::hypot(xi, sigma) };
                    const double wKModel{ std::fma(b, rho * xi + R, a) };

                    // Calculate SSE between model and market variance
                    const double r{ wKModel - wKMarket };
                    SSE += r * r;

                    if (grad)
                    {
                        // Precompute variables
                        const double invR{ 1.0 / R };

                        // Calculate partial derivatives
                        const double dwda{ 1.0 };                     // ∂w/∂a
                        const double dwdb{ (rho * xi + R) };          // ∂w/∂b
                        const double dwdrho{ (b * xi) };              // ∂w/∂ρ
                        const double dwdm{ -b * (rho + xi * invR) };  // ∂w/∂m
                        const double dwdsigma{ b * (sigma * invR) };  // ∂w/∂σ

                        // ∇f += 2 * r * ∂w/∂θ  (variance-SSE)
                        g[0] += 2.0 * r * dwda;
                        g[1] += 2.0 * r * dwdb;
                        g[2] += 2.0 * r * dwdrho;
                        g[3] += 2.0 * r * dwdm;
                        g[4] += 2.0 * r * dwdsigma;
                    }
                }

                if (grad)
                {
                    // Provide gradient to NLopt
                    const auto mCopy = (std::min)(static_cast<std::size_t>(n), g.size());
                    std::copy_n(g.data(), mCopy, grad);
                }
                return SSE;
            },
            const_cast<ObjCtx*>(&obj) // nlopt needs void*
        );
    }

    template <::nlopt::algorithm Algo>
    void evalCal(const Params& sviSlice,
        const opt::nlopt::Optimizer<5, Algo>& optimizer,
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

        // ---- Calendar no-arb: w_curr(k) >= w_prev(k) + eps ----
        std::size_t nCalViol = 0;
        double minMargin = std::numeric_limits<double>::infinity();
        double kAtWorst = 0.0;

        for (std::size_t i = 0; i < kSlice.size(); ++i)
        {
            const double k = kSlice[i];
            const double wPrev = wKPrevSlice[i] + eps;
            const double wCurr = wk(a, b, rho, m, sigma, k);
            const double margin = wCurr - wPrev; // should be >= 0
            if (margin < minMargin) { minMargin = margin; kAtWorst = k; }
            if (margin < -tol) ++nCalViol;
        }

        UV_WARN(minMargin < -tol,
            std::format("No-arbitrage: calendar violated: "
                "min[w_curr - (w_prev+eps)] = {:.6e} at k = {:.4f} "
                "(violations {} / {}, tol = {:.2e}, eps = {:.2e})",
                minMargin, kAtWorst, nCalViol, kSlice.size(),
                tol, eps));
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
