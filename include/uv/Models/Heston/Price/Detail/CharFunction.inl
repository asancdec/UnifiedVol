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

#pragma once

#include "Base/Types.hpp"

#include <cmath>
#include <concepts>
#include <limits>

namespace uv::models::heston::price::detail
{

template <std::floating_point T> Complex<T> charFunction(
    T kappa,
    T kappaThetaDivSigma2,
    T sigma2,
    T v0,
    T t,
    Complex<T> tDivTwo,
    Complex<T> sigmaRho,
    Complex<T> u
) noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};

    const Complex<T> beta{kappa + sigmaRho * u};

    const Complex<T> uu{u * (u + i)};

    const Complex<T> sigma2uu{sigma2 * uu};

    const Complex<T> D{std::sqrt(beta * beta + sigma2uu)};

    const Complex<T> r{
        (std::real(beta * std::conj(D)) > T{0}) ? -sigma2uu / (beta + D) : beta - D
    };

    const Complex<T> DT{D * t};

    Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T{1} + std::abs(DT)))
            ? math::expm1Complex(-DT) * T{0.5} / D
            : tDivTwo
    };

    const Complex<T> ry{-r * y};

    return kappaThetaDivSigma2 * (r * t - T{2} * math::log1pComplex<T>(ry)) +
           v0 * (uu * y / (T{1} + ry));
}

template <std::floating_point T> [[gnu::hot]] CharFunCache<T> charFunctionCached(
    T kappa,
    T kappaThetaDivSigma2,
    T sigma2,
    T v0,
    T t,
    Complex<T> tDivTwo,
    Complex<T> sigmaRho,
    Complex<T> u
) noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};

    const Complex<T> beta{kappa + sigmaRho * u};

    const Complex<T> uu{u * (u + i)};

    const Complex<T> sigma2uu{sigma2 * uu};

    const Complex<T> D{std::sqrt(beta * beta + sigma2uu)};

    const Complex<T> betaPlusD{beta + D};

    const Complex<T> betaMinusD{
        (std::real(beta * std::conj(D)) > T{0}) ? -sigma2uu / betaPlusD : beta - D
    };

    const Complex<T> DT{D * t};

    const Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T{1} + std::abs(DT)))
            ? math::expm1Complex(-DT) * T{0.5} / D
            : tDivTwo
    };

    const Complex<T> ry{-betaMinusD * y};

    const Complex<T> A{
        kappaThetaDivSigma2 * (betaMinusD * t - T{2} * math::log1pComplex<T>(ry))
    };

    const Complex<T> B{uu * y / (T{1} + ry)};

    const Complex<T> eDT{std::exp(-DT)};

    const Complex<T> g{betaMinusD / betaPlusD};

    const Complex<T> Q{T{1} - g * eDT};

    const Complex<T> invQ{T{1} / Q};

    const Complex<T> R{T{1} - g};

    return CharFunCache{
        A + v0 * B,
        A,
        B,
        beta,
        D,
        betaMinusD,
        uu,
        eDT,
        g,
        Q,
        invQ,
        R,
        betaMinusD * t - T{2} * std::log(Q / R),
        (T{1} - eDT) * invQ,
        betaPlusD * betaPlusD
    };
}

} // namespace uv::models::heston::price::detail

#include "Models/Heston/Price/Detail/CharFunction.inl"