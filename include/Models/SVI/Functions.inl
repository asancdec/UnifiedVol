// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Inline implementations of performance-critical SVI routines.
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

#include "Core/Functions.hpp"
#include "Core/Matrix/Functions.hpp"
#include "Math/Interpolation/Interpolator.hpp"
#include "Utils/Aux/Errors.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace uv::models::svi
{

template <std::floating_point T, ::nlopt::algorithm Algo>
Vector<Params<T>> calibrate(
    const core::VolSurface<T>& volSurface,
    const opt::nlopt::Optimizer<4, Algo>& prototype
)
{
    return calibrate<Real, Algo>(
        volSurface.tenors(),
        volSurface.logKFMatrix(),
        volSurface.totVarMatrix(),
        prototype
    );
}

template <std::floating_point T, ::nlopt::algorithm Algo>
Vector<Params<T>> calibrate(
    const Vector<T>& tenors,
    const core::Matrix<T>& kMatrix,
    const core::Matrix<T>& wMMatrix,
    const opt::nlopt::Optimizer<4, Algo>& prototype
)
{
    // ---------- Validate inputs ----------

    detail::validateInputs<T>(tenors, kMatrix, wMMatrix);

    // ---------- Extract data ----------

    // NLopt API accepts double type only
    const core::Matrix<double> kMatrixD{core::convertMatrix<double>(kMatrix)};
    const core::Matrix<double> wMMatrixD{core::convertMatrix<double>(wMMatrix)};

    const std::size_t numTenors{tenors.size()};
    const std::size_t numStrikes{kMatrixD.cols()};

    // ---------- Initialize params container ----------

    Vector<Params<T>> params;
    params.reserve(numTenors);

    // ---------- Calibrate per slice ----------

    for (std::size_t i = 0; i < numTenors; ++i)
    {
        const Params<T>* prev = (i == 0) ? nullptr : &params.back();

        params.emplace_back(detail::calibrateSlice<T>(
            tenors[i],
            kMatrixD[i],
            wMMatrixD[i],
            prototype,
            prev,
            numStrikes
        ));
    }

    return params;
}

template <std::floating_point T>
core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Vector<Params<T>>& params)
{
    // ---------- Extract data ----------

    const std::size_t numTenors{volSurface.numTenors()};
    const std::size_t numStrikes{volSurface.numStrikes()};
    const core::Matrix<T>& kMatrix{volSurface.logKFMatrix()};

    // ---------- Validate inputs ----------

    UV_REQUIRE(
        params.size() == numTenors,
        ErrorCode::InvalidArgument,
        "buildSurface: params size must equal number of tenors"
    );

    // ---------- Build surface  ----------

    // Copy surface
    core::VolSurface<T> sviVolSurface{volSurface};

    // Set total variance
    sviVolSurface.setTotVar(core::generateIndexed<Real>(
        numTenors,
        numStrikes,
        [&](std::size_t i, std::size_t j)
        {
            std::span<const T> kMatrixRow{kMatrix[i]};
            const Params<T>& p{params[i]};

            return detail::calculateWk(p.a, p.b, p.rho, p.m, p.sigma, kMatrixRow[j]);
        }
    ));

    return sviVolSurface;
}

template <std::floating_point T> T gk(T a, T b, T rho, T m, T sigma, T k) noexcept
{
    T x{k - m};                    // x := k-m
    T sigmaSquared{sigma * sigma}; // sigma^2

    // R := sqrt(x^2 + sigma^2)
    T R{std::sqrt(std::fma(x, x, sigmaSquared))};
    T invR{1.0 / R}; // invR := 1 / R

    T wk{std::fma(b, (rho * x + R), a)}; // w(k) = a + b*(rho*x + R)
    T wkInv{1.0 / wk};                   // 1 / w(k)

    T wkD1{b * (rho + x * invR)}; // w'(k) = b * (rho + x/R)
    T wkD1Squared{wkD1 * wkD1};   // w'(k)^2

    T invR2{invR * invR};      // 1 / R^2
    T invRCubed{invR2 * invR}; // 1 / R^3

    T wkD2{b * sigmaSquared * invRCubed}; // w''(k) = b * sigma^2 / R^3

    T A{1.0 - 0.5 * k * wkD1 * wkInv}; // A := 1 - k * w'/(2 * w)
    T B{wkInv + 0.25};                 // B := 1/w(k) + 1/4

    return (A * A) - 0.25 * wkD1Squared * B + wkD2 * 0.5;
}

} // namespace uv::models::svi

namespace uv::models::svi::detail
{
struct ObjectiveContexts
{
    const double* k;    // Pointer to log-forward moneyness
    const double* wM;   // Pointer to market total variance
    std::size_t n;      // Number of points
    const double atmWK; // Atm variance
};

struct CalendarContexts
{
    double k;      // Log-forward moneyness
    double prevWk; // Total variance of the previous slice
    double eps;    // Epsilon value
    double atmWK;  // Atm variance
};

struct ConvexityContexts
{
    double k;     // Log-forward moneyness
    double atmWK; // Atm variance
};

struct GkCache
{
    double x;            // x := k-m
    double R;            // R:= sqrt(x^2 + sigma^2)
    double invR;         // invR := 1 / R
    double wk;           // w(k) = a + b*(rho*x + R)
    double wkD1;         // w'(k) = b * (rho + x/R)
    double wkD1Squared;  // w'(k)^2
    double invRCubed;    // 1/(R^3)
    double invR5;        // 1/(R^5)
    double sigmaSquared; // sigma^2
    double wkD2;         // w''(k) = b * sigma^2 / R^3
    double A;            // A := 1 - k * w'/(2 * w)
    double B;            // B := 1/w(k) + 1/4
    double wkInv;        // 1/w(k)
    double wkSquaredInv; // 1/w(k)^2
    double R0;           // R0 := sqrt(m^2 + sigma^2)
    double invR0;        // 1/R0

    explicit GkCache(
        double a,
        double b,
        double rho,
        double m,
        double sigma,
        double k
    ) noexcept
    {
        sigmaSquared = sigma * sigma; // sigma^2

        R0 = std::sqrt(std::fma(m, m, sigmaSquared)); // R0 := sqrt(m^2 + sigma^2)
        invR0 = 1.0 / R0;                             // 1/R0

        x = k - m;                                   // x := k-m
        R = std::sqrt(std::fma(x, x, sigmaSquared)); // R:= sqrt(x^2 + sigma^2)
        invR = 1.0 / R;                              // invR := 1 / R

        wk = std::fma(b, (rho * x + R), a); // w(k) = a + b*(rho*x + R)
        wkInv = 1.0 / wk;                   // 1/w(k)
        wkSquaredInv = wkInv * wkInv;       // 1/w(k)^2

        wkD1 = b * (rho + x / R);  // w'(k) = b * (rho + x/R)
        wkD1Squared = wkD1 * wkD1; // w'(k)^2

        invRCubed = invR * invR * invR;  // 1/(R^3)
        invR5 = invRCubed * invR * invR; // 1/(R^5)

        wkD2 = b * sigmaSquared * invRCubed; // w''(k) = b * sigma^2 / R^3

        A = 1.0 - 0.5 * k * wkD1 * wkInv; // A := 1 - k * w'/(2 * w)
        B = wkInv + 0.25;                 // B := 1/w(k) + 1/4
    }
};

template <std::floating_point T>
void validateInputs(
    const Vector<T>& tenors,
    const core::Matrix<T>& kMatrix,
    const core::Matrix<T>& wMMatrix
)
{
    UV_REQUIRE(
        !tenors.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: tenors is empty"
    );

    UV_REQUIRE(
        !kMatrix.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: kMatrix is empty"
    );

    UV_REQUIRE(
        !wMMatrix.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: wMMatrix is empty"
    );

    UV_REQUIRE(
        kMatrix.rows() == tenors.size(),
        ErrorCode::InvalidArgument,
        "validateInputs: kMatrix rows must equal number of tenors"
    );

    UV_REQUIRE(
        wMMatrix.rows() == tenors.size(),
        ErrorCode::InvalidArgument,
        "validateInputs: wMMatrix rows must equal number of tenors"
    );

    const std::size_t numTenors{tenors.size()};
    const std::size_t numStrikes{kMatrix.cols()};

    UV_REQUIRE(
        wMMatrix.cols() == numStrikes,
        ErrorCode::InvalidArgument,
        "validateInputs: wMMatrix columns must equal kMatrix columns"
    );

    for (std::size_t i = 1; i < numTenors; ++i)
    {
        UV_REQUIRE(
            tenors[i] > tenors[i - 1],
            ErrorCode::InvalidArgument,
            "validateInputs: tenors must be strictly increasing"
        );
    }

    for (std::size_t i = 0; i < numTenors; ++i)
    {
        for (std::size_t j = 1; j < numStrikes; ++j)
        {
            UV_REQUIRE(
                kMatrix[i][j] > kMatrix[i][j - 1],
                ErrorCode::InvalidArgument,
                "validateInputs: kMatrix rows must be strictly increasing"
            );
        }
    }
}

template <std::floating_point T, ::nlopt::algorithm Algo>
Params<T> calibrateSlice(
    const T t,
    std::span<const double> kSlice,
    std::span<const double> wKSlice,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    const Params<T>* prevParams,
    const std::size_t numStrikes
)
{
    // ---------- Calculate  ----------

    const double atmWK{math::interp::PchipInterpolator<double>{}(0.0, kSlice, wKSlice)};
    const double logKFMin{core::minValue(kSlice)};
    const double logKFMax{core::maxValue(kSlice)};

    // ---------- Clone  ----------

    opt::nlopt::Optimizer optimizer{prototype.fresh()};
    optimizer.setUserValue(atmWK);

    // ---------- Allocate  ----------

    // NOTE: Must outlive optimizer.optimize();
    // NLopt holds pointers to elements.
    Vector<CalendarContexts> calendarContexts{};

    // ---------- Configure ----------

    if (!prevParams)
    {
        // Bounds and guess
        optimizer
            .setGuessBounds(coldGuess(), lowerBounds(logKFMin), upperBounds(logKFMax));
    }
    else
    {
        // Bounds and guess
        optimizer.setGuessBounds(
            warmGuess(*prevParams),
            lowerBounds(logKFMin),
            upperBounds(logKFMax)
        );

        // Calendar constraints
        calendarContexts.resize(numStrikes + 2);

        fillCalendarContexts<T, Algo>(
            calendarContexts,
            numStrikes,
            *prevParams,
            optimizer,
            kSlice,
            atmWK,
            logKFMin,
            logKFMax
        );

        addCalendarConstraint(optimizer, calendarContexts);
    }

    // ---------- Standard constraints ----------

    addWMinConstraint(optimizer);
    addMinSlopeConstraint(optimizer);
    addMaxSlopeConstraint(optimizer);

    // ---------- Convexity constraints ----------

    Vector<ConvexityContexts> convexityCtxs;
    convexityCtxs.reserve(numStrikes);
    for (double k : kSlice)
    {
        convexityCtxs.push_back({k, atmWK});
        addConvexityConstraint(optimizer, convexityCtxs.back());
    }

    // ---------- Set objective function ----------

    ObjectiveContexts obj{
        kSlice.data(),
        wKSlice.data(),
        numStrikes,
        optimizer.userValue()
    };

    setMinObjective(optimizer, obj);

    // ---------- Run optimization ----------

    Vector<double> params{optimizer.optimize()};

    return Params<T>{
        t,
        aParam(atmWK, T(params[0]), T(params[1]), T(params[2]), T(params[3])),
        T(params[0]),
        T(params[1]),
        T(params[2]),
        T(params[3])
    };
}

template <std::floating_point T>
std::array<double, 4> warmGuess(const Params<T>& params) noexcept
{
    return {
        // Explicit conversion
        double(params.b),    // b
        double(params.rho),  // rho
        double(params.m),    // m
        double(params.sigma) // sigma
    };
}

template <::nlopt::algorithm Algo>
void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{
    optimizer.addInequalityConstraint(
        +[](unsigned /*n*/, const double* x, double* grad, void* data) noexcept -> double
        {
            const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

            // Extract data
            const double eps{opt.eps()};
            const double atmWK{opt.userValue()};
            const double b{x[0]};
            const double rho{x[1]};
            const double m{x[2]};
            const double sigma{x[3]};

            // Precomputes
            const double s2{sigma * sigma};
            const double R{std::sqrt(m * m + s2)};
            const double S{std::sqrt(1.0 - rho * rho)};

            // ATM base and a
            const double base{-rho * m + R};
            const double a{atmWK - b * base};

            // wMin = a + b * sigma * S
            const double wMin{std::fma(b, sigma * S, a)};

            if (grad)
            {
                // same formulas you had
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
void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{
    // ---------- Right wing: b*(1 + rho) >= eps ----------

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
        {
            // Reinterpret opaque pointer
            const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

            // Extract data
            const double eps{opt.eps()};
            const double b{x[0]};
            const double onePlusRho{1.0 + x[1]};

            if (grad)
            {
                grad[0] = -onePlusRho;
                grad[1] = -b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return eps - b * onePlusRho;
        },
        &optimizer
    );

    // ---------- Left wing: b*(1 - rho) >= eps ----------

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
        {
            // Reinterpret opaque pointer
            const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

            // Extract data
            const double eps{opt.eps()};
            const double b{x[0]};
            const double oneMinusRho{1.0 - x[1]};

            if (grad)
            {
                grad[0] = -oneMinusRho;
                grad[1] = b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return eps - b * oneMinusRho;
        },
        &optimizer
    );
}

template <std::floating_point T>
CalendarContexts
genCalendarContext(double k, const Params<T>& params, double eps, double atmWK) noexcept
{
    return CalendarContexts{
        k,
        calculateWk(
            // Avoid implicit conversions
            double(params.a),
            double(params.b),
            double(params.rho),
            double(params.m),
            double(params.sigma),
            k
        ),
        eps,
        atmWK
    };
}

template <std::floating_point T, ::nlopt::algorithm Algo>
void fillCalendarContexts(
    Vector<CalendarContexts>& calendarContexts,
    std::size_t numStrikes,
    const Params<T>& params,
    const opt::nlopt::Optimizer<4, Algo>& optimizer,
    std::span<const double> kSlice,
    double atmWK,
    double logKFMin,
    double logKFMax,
    double delta
) noexcept
{
    // ---------- Extract  ----------

    const double eps{optimizer.eps()};

    // ---------- Fill  ----------

    // Market strikes
    for (std::size_t i = 0; i < numStrikes; ++i)
    {
        calendarContexts[i] = (genCalendarContext<T>(kSlice[i], params, eps, atmWK));
    }

    // Strikes outside market range

    // Lower strike
    calendarContexts[numStrikes] =
        (genCalendarContext<T>(logKFMin - delta, params, eps, atmWK));

    // Upper strike
    calendarContexts[numStrikes + 1] =
        (genCalendarContext<T>(logKFMax + delta, params, eps, atmWK));
}

template <::nlopt::algorithm Algo>
void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{
    // ---------- Right wing: b*(1 + rho) <= 2 ----------

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) noexcept -> double
        {
            // Extract data
            const double b{x[0]};
            const double rho{x[1]};

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
        +[](unsigned, const double* x, double* grad, void*) noexcept -> double
        {
            // Extract data
            const double b{x[0]};
            const double rho{x[1]};

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
void addCalendarConstraint(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    Vector<CalendarContexts>& ctx
) noexcept
{
    for (auto& ctxK : ctx)
    {
        optimizer.addInequalityConstraint(
            +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
            {
                const CalendarContexts& ctxK{*static_cast<const CalendarContexts*>(data)};

                // Extract data
                const double atmWK{ctxK.atmWK};
                const double k{ctxK.k};
                const double prevWk{ctxK.prevWk};
                const double eps{ctxK.eps};
                const double b{x[0]};
                const double rho{x[1]};
                const double m{x[2]};
                const double sigma{x[3]};

                // Precomputes
                const double xi{k - m};
                const double s2{sigma * sigma};

                // ATM pieces
                const double R0{std::sqrt(m * m + s2)};
                const double base{-rho * m + R0};
                const double a{atmWK - b * base};

                // Strike-specific
                const double Rk{std::sqrt(xi * xi + s2)};
                const double wK{std::fma(b, (rho * xi + Rk), a)};

                // ----------- No-grad path -----------

                if (!grad)
                    return prevWk + eps - wK;

                // ----------- Grad path -----------
                const double invRk{1.0 / Rk};
                const double invR0{1.0 / R0};
                const double mInvR0{m * invR0};
                const double sigmaInvR0{sigma * invR0};

                grad[0] = -((rho * xi + Rk) - base);
                grad[1] = -(b * (xi + m));
                grad[2] = (b * (mInvR0 + xi * invRk));
                grad[3] = -(b * (sigma * invRk - sigmaInvR0));

                return prevWk + eps - wK;
            },
            &ctxK
        );
    }
}

template <::nlopt::algorithm Algo>
void addConvexityConstraint(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    ConvexityContexts& ctx
) noexcept
{
    optimizer.addInequalityConstraint(
        +[](unsigned /*n*/, const double* x, double* grad, void* data) noexcept -> double
        {
            const auto& ctx = *static_cast<const ConvexityContexts*>(data);

            // Extract data
            const double k{ctx.k};
            const double atmWK{ctx.atmWK};
            const double b{x[0]};
            const double rho{x[1]};
            const double m{x[2]};
            const double sigma{x[3]};

            // ATM precomputes
            const double s2{sigma * sigma};
            const double R0{std::sqrt(m * m + s2)};
            const double base{-rho * m + R0};
            const double a{atmWK - b * base};

            // Cache once
            const GkCache p{a, b, rho, m, sigma, k};
            const double g{gk(p)};

            // ----------- No-grad path -----------

            if (!grad)
                return -g;

            // ----------- Grad path -----------

            const std::array<double, 4> dg{gkGrad(b, rho, m, sigma, k, p)};

            grad[0] = -dg[0];
            grad[1] = -dg[1];
            grad[2] = -dg[2];
            grad[3] = -dg[3];

            return -g;
        },
        &ctx
    );
}

template <::nlopt::algorithm Algo>
void setMinObjective(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const ObjectiveContexts& ctx
) noexcept
{
    optimizer.setMinObjective(
        +[](unsigned /*n*/, const double* x, double* grad, void* data) noexcept -> double
        {
            const ObjectiveContexts& ctx = *static_cast<const ObjectiveContexts*>(data);

            // Extract data
            const double atmWK{ctx.atmWK};
            const std::size_t N{ctx.n};
            const double b{x[0]};
            const double rho{x[1]};
            const double m{x[2]};
            const double sigma{x[3]};

            // Precomputes (ATM)
            const double s2{sigma * sigma};
            const double R0{std::sqrt(m * m + s2)};
            const double invR0{1.0 / R0};
            const double base{-rho * m + R0};
            const double a{atmWK - b * base};

            // Initialize
            double SSE{0.0};

            // ----------- No-grad path -----------
            if (!grad)
            {
                for (std::size_t i = 0; i < N; ++i)
                {
                    const double k{ctx.k[i]};
                    const double wM{ctx.wM[i]};

                    const double xi{k - m};
                    const double R{std::sqrt(xi * xi + s2)};

                    const double wK{std::fma(b, (rho * xi + R), a)};
                    const double r{wK - wM};

                    SSE = std::fma(r, r, SSE);
                }
                return SSE;
            }

            // ----------- Grad path -----------

            double g0{0.0}, g1{0.0}, g2{0.0}, g3{0.0};
            const double sigmaInvR0{sigma * invR0};
            const double mInvR0{m * invR0};

            for (std::size_t i = 0; i < N; ++i)
            {
                const double k{ctx.k[i]};
                const double wM{ctx.wM[i]};

                const double xi{k - m};
                const double R{std::sqrt(xi * xi + s2)};
                const double invR{1.0 / R};

                const double wK{std::fma(b, (rho * xi + R), a)};
                const double r{wK - wM};

                SSE = std::fma(r, r, SSE);

                // keep your exact algebra
                g0 += 2.0 * r * ((rho * xi + R) - base);
                g1 += 2.0 * r * (b * (xi + m));
                g2 += 2.0 * r * (-b * (mInvR0 + xi * invR));
                g3 += 2.0 * r * (b * (sigma * invR - sigmaInvR0));
            }

            grad[0] = g0;
            grad[1] = g1;
            grad[2] = g2;
            grad[3] = g3;

            return SSE;
        },
        const_cast<ObjectiveContexts*>(&ctx)
    );
}

template <std::floating_point T>
T calculateWk(T a, T b, T rho, T m, T sigma, T k) noexcept
{
    const T x{k - m};
    const T R{std::sqrt(std::fma(x, x, sigma * sigma))};
    return std::fma(b, (rho * x + R), a);
}
} // namespace uv::models::svi::detail