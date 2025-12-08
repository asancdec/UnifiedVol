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

#include <cmath>
#include <numbers>

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
        constexpr T invSqrt2Pi{ T(1.0) / std::sqrt(T(2.0) * std::numbers::pi_v<T>) };
        return invSqrt2Pi * std::exp(-T(0.5) * x * x);
    }

    template <std::floating_point T>
    T d1BS(T t,
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
        T d1{ d1BS(t, r, q, vol, S, K)};
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
}  // namespace uv::math