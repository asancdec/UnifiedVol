// SPDX-License-Identifier: Apache-2.0
/*
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

#include <cmath>
#include <complex>
#include <numbers>

namespace uv::math
{
template <std::floating_point T> Complex<T> log1pComplex(Complex<T> z) noexcept
{

    const T a{std::real(z)};
    const T b{std::imag(z)};

    if (std::abs(a) < 0.5 && std::abs(b) < 0.5)
    {

        return {
            0.5 * std::log1p(std::fma(a, a, std::fma(2.0, a, b * b))),
            std::atan2(b, 1.0 + a)
        };
    }

    return std::log(Complex<T>(1.0 + a, b));
}

template <std::floating_point T> T cosm1(T b) noexcept
{

    const T s{std::sin(b * T{0.5})};
    return -2.0 * s * s;
}

template <std::floating_point T> Complex<T> expm1Complex(Complex<T> z) noexcept
{

    if (std::abs(z) < 1.0)
    {

        const T a{std::real(z)};
        const T b{std::imag(z)};

        const T cm1{cosm1(b)};

        const T em1{std::expm1(a)};

        return {em1 * (cm1 + 1.0) + cm1, std::sin(b) * std::exp(a)};
    }

    return std::exp(z) - T(1);
}

template <std::floating_point T> T normalCDF(T x) noexcept
{
    return std::erfc(-x / std::sqrt(2.0)) * 0.5;
}

template <std::floating_point T> T normalPDF(T x) noexcept
{
    constexpr T invSqrt2Pi = std::numbers::inv_sqrtpi_v<T> / std::numbers::sqrt2_v<T>;
    return invSqrt2Pi * std::exp(-0.5 * x * x);
}
} // namespace uv::math