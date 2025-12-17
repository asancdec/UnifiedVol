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
#include "Utils/Aux/Errors.hpp"

#include <cmath>
#include <numbers>
#include <string>
#include <format>
#include <algorithm>

namespace uv::math
{
    template <std::floating_point T>
    Complex<T> log1pComplex(const Complex<T>& z) noexcept
    {
        // a := Re(z),  b := Im(z)
        const T a{ std::real(z) };
        const T b{ std::imag(z) };

        // Small-box branch: |a|, |b| < 1/2
        if (std::abs(a) < T(0.5) && std::abs(b) < T(0.5))
        {
            // Re := 0.5 * log1p(a^2 + 2a + b^2)
            // Im := atan2(b, 1 + a)
            return
            {
                T(0.5) * std::log1p(std::fma(a, a, std::fma(T(2), a, b * b))),
                std::atan2(b, T(1) + a)
            };
        }

        // General case: ln(1 + z)
        return std::log(Complex<T>(T(1) + a, b));
    }

    template <std::floating_point T>
    T cosm1(T b) noexcept 
    {
        // cosm1(b) := cos(b) - 1 = -2 sin²(b / 2)
        const T s{ std::sin(b * T(0.5)) };
        return T(-2) * s * s;                
    }

    template <std::floating_point T>
    Complex<T> expm1Complex(const Complex<T>& z) noexcept 
    {   
        // If |z| < 1 : use numerically stable expansion
        if (std::abs(z) < T(1))
        {
            // a := Re(z),  b := Im(z)
            const T a{ std::real(z) };
            const T b{ std::imag(z) };

            // cosm1(b) := cos(b) - 1 (accurate for small b)
            const T cm1{ cosm1(b) };

            // em1 := e^a - 1 (accurate near 0)
            const T em1{ std::expm1(a) };

            // Re := em1 * (cm1 + 1) + cm1
            // Im := sin(b) * e^a
            return { em1 * (cm1 + T(1)) + cm1, std::sin(b) * std::exp(a) };
        }

        // Otherwise use direct definition e^z - 1
        return std::exp(z) - T(1);
    }

    template <std::floating_point T>
    T normalCDF(T x) noexcept
    {
        return std::erfc(-x / std::sqrt(T(2.0))) * T(0.5);
    }

    template <std::floating_point T>
    T normalPDF(T x) noexcept
    {
        constexpr T invSqrt2Pi = std::numbers::inv_sqrtpi_v<T> / std::numbers::sqrt2_v<T>;
        return invSqrt2Pi * std::exp(-T(0.5) * x * x);
    }

    template <std::floating_point T>
    T d1BlackScholes(T t,
        T r,
        T q,
        T vol,
        T S,
        T K) noexcept
    {
        return std::fma(t, (r - q + vol * vol * T(0.5)), std::log(S / K)) / (std::sqrt(t) * vol);
    }

    template <std::floating_point T>
    T blackScholes(T t,
        T r,
        T q,
        T vol,
        T S,
        T K,
        bool isCall) noexcept
    {
        T d1{ d1BlackScholes(t, r, q, vol, S, K)};
        T d2{ std::fma(-vol, std::sqrt(t), d1) };

        if (isCall)
        {
            return S * std::exp(-q * t) * normalCDF(d1) - K * std::exp(-r * t) * normalCDF(d2);
        }
        else
        {
            return K * std::exp(-r * t) * normalCDF(-d2) - S * std::exp(-q * t) * normalCDF(-d1);
        }
    }

    template <std::floating_point T>
    T vegaBS(T d1,
        T t,
        T q,
        T S) noexcept
    {
        return S * std::exp(-q * t) * normalPDF(d1) * std::sqrt(t);
    }

    template <std::floating_point T>
    T volgaBS(T vega,
        T d1,
        T t,
        T vol) noexcept
    {
        T d2{ std::fma(-vol, std::sqrt(t), d1) };
        return vega * (d1 * d2) / vol;
    }

    template <std::floating_point T>
    T impliedVolBS(T callPrice,
        T t,
        T r,
        T q,
        T S,
        T K)
    {
        // ---------- Optimization parameters ----------

        const T TOL{ 1e-14 };
        const unsigned int EVAL{ 100 };

        // ---------- Sanity checks ----------

        // Spot must be positive
        UV_REQUIRE(S > T(0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: spot S = {:.6f} must be > 0", S)
        );

        // Strike must be positive
        UV_REQUIRE(K > T(0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: strike K = {:.6f} must be > 0", K)
        );

        // Tenor must be positive
        UV_REQUIRE(t > T(0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: tenor t = {:.6f} must be > 0", t)
        );

        // Market price must be non-negative
        UV_REQUIRE(callPrice >= T(0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: market price = {:.6f} must be >= 0", callPrice)
        );

        // ---------- Initial guess ----------

        // Log(F/K)
        const T logKF{ std::log((S * std::exp((r - q) * t)) / K) };

        // Heuristic guess
        const T volGuess
        {
            (std::fabs(logKF) < T(1e-6))
                ? T(0.3)
                : std::sqrt(T(2) * std::fabs(logKF) / t)
        };

        // Clamp bounds
        const T volMin{ T(1e-4) };
        const T volMax{ T(5.0) };

        // Clamp
        T vol{ std::clamp(volGuess, volMin, volMax) };

        // Warn if clamped
        UV_WARN(vol != volGuess,
            std::format(
                "impliedVolBS: initial guess = {:.6f} clamped to {:.6f} (lb = {:.4f}, ub = {:.4f})",
                volGuess, vol, volMin, volMax
            )
        );

        // ---------- Halley's method ----------

        for (unsigned int i = 0U; i < EVAL; ++i)
        {

            const T priceBS{ blackScholes(t, r, q, vol, S, K) };
            const T objEval{ priceBS - callPrice };

            // Check Absolute and relative tolerance threshold
            if (std::fabs(objEval) < TOL * (T(1) + priceBS)) break;

            // Derivatives
            const T d1{ d1BlackScholes<T>(t, r, q, vol, S, K) };
            const T vega{ vegaBS<T>(d1, t, q, S) };
            const T volga{ volgaBS<T>(vega, d1, t, vol) };

            // Update volatility using Halley's method
            vol -= (T(2) * objEval * vega) / (T(2) * (vega * vega) - objEval * volga);
        }

        // ---------- Evaluate calibration ----------

        // Throw if volatility is negative or larger than 100
        UV_REQUIRE(
            (vol > T(0)) && (vol < T(100)),
            ErrorCode::CalibrationError,
            std::format(
                "impliedVolBS: resulting volatility out of bounds: vol = {:.6f} "
                "(lb = {:.4f}, ub = {:.4f})",
                vol, T(0.0), T(100)
            )
        );

        const T finalPrice{ blackScholes<T>(t, r, q, vol, S, K) };
        const T finalResidual{ finalPrice - callPrice };
        const T finalTol{ TOL * (T(1) + std::fabs(finalPrice)) };

        // Warn if no convergence
        UV_WARN(std::fabs(finalResidual) > finalTol,
            std::format(
                "impliedVolBS: no convergence after {} iterations "
                "(|f| = {:.3e} > tol = {:.3e}, vol = {:.6f})",
                EVAL, std::fabs(finalResidual), finalTol, vol
            )
        );

        return vol;
    }

    template <std::floating_point T>
    Vector<T> impliedVolBS(const Vector<T>& callPrices,
        T t,
        T r,
        T q,
        T S,
        const Vector<T>& K
    )
    {
        // ---------- Validate inputs ----------

        const std::size_t N{ callPrices.size() };

        UV_REQUIRE(
            K.size() == N,
            ErrorCode::InvalidArgument,
            "impliedVolBS: price / strike size mismatch"
        );

        // ---------- Compute implied volatilities ----------

        Vector<T> out(N);

        for (std::size_t i = 0; i < N; ++i)
        {
            out[i] = math::impliedVolBS(
                callPrices[i],
                t,
                r,
                q,
                S,
                K[i]
            );
        }

        return out;
    }

    template <std::floating_point T>
    core::Matrix<T> impliedVolBS(const core::Matrix<T>& callPrices,
        const Vector<T>& t,
        const Vector<T>& r,
        const Vector<T>& q,
        T S,
        const Vector<T>& K
    )
    {
        // ---------- Validate inputs ----------

        const std::size_t Nt{ t.size() };
        const std::size_t Nk{ K.size() };

        UV_REQUIRE(
            r.size() == Nt &&
            q.size() == Nt &&
            callPrices.rows() == Nt &&
            callPrices.cols() == Nk,
            ErrorCode::InvalidArgument,
            "impliedVolBS(matrix): input size mismatch"
        );

        // ---------- Compute implied volatilities ----------

        core::Matrix<T> out(Nt, Nk);

        for (std::size_t i = 0; i < Nt; ++i)
        {
            std::span<T> outRow{out[i]};
            std::span<const T> callPricesRow{callPrices[i]};
            const T ti{ t[i] };
            const T ri{ r[i] };
            const T qi{ q[i] };

            for (std::size_t j = 0; j < Nk; ++j)
            {
                outRow[j] = math::impliedVolBS(
                    callPricesRow[j],
                    ti,
                    ri,
                    qi,
                    S,
                    K[j]
                );
            }
        }

        return out;
    }

}  // namespace uv::math