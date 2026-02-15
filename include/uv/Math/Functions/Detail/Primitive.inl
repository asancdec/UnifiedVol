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

template <std::floating_point T> Complex<T> invComplex(Complex<T> z) noexcept
{

    const T invNorm{T{1} / std::norm(z)};
    return std::conj(z) * invNorm;
}

template <std::floating_point T> Complex<T> log1pComplex(Complex<T> z) noexcept
{

    const T a{std::real(z)};
    const T b{std::imag(z)};

    const T ap1{T{1} + a};

    const T aTan2ap1{std::atan2(b, ap1)};

    if (std::abs(a) < T{0.5} && std::abs(b) < T{0.5})
    {

        return {T{0.5} * std::log1p(a * (a + T{2}) + b * b), aTan2ap1};
    }

    return {std::log(std::hypot(ap1, b)), aTan2ap1};
}

template <std::floating_point T> T cosm1(T b) noexcept
{

    const T s{std::sin(b * T{0.5})};
    return T{-2} * s * s;
}

template <std::floating_point T> Complex<T> expm1Complex(Complex<T> z) noexcept
{

    if (std::norm(z) < T{1})
    {

        const T a{std::real(z)};
        const T b{std::imag(z)};

        const T cm1{cosm1(b)};

        const T em1{std::expm1(a)};

        return {em1 * (cm1 + T{1}) + cm1, std::sin(b) * (em1 + T{1})};
    }

    return std::exp(z) - T(1);
}

template <std::floating_point T> T normalCDF(T x) noexcept
{
    return std::erfc(-x / std::sqrt(T{2})) * T{0.5};
}

template <std::floating_point T> T normalPDF(T x) noexcept
{
    constexpr T invSqrt2Pi = std::numbers::inv_sqrtpi_v<T> / std::numbers::sqrt2_v<T>;
    return invSqrt2Pi * std::exp(-T{0.5} * x * x);
}
} // namespace uv::math