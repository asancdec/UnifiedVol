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
#include "Math/Functions/Primitive.hpp"

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
        (std::real(beta * std::conj(D)) > T{0}) ? -sigma2uu * math::invComplex(beta + D)
                                                : beta - D
    };

    const Complex<T> DT{D * t};

    Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T{1} + std::abs(DT)))
            ? math::expm1Complex(-DT) * T{0.5} * math::invComplex(D)
            : tDivTwo
    };

    const Complex<T> ry{-r * y};

    return kappaThetaDivSigma2 * (r * t - T{2} * math::log1pComplex<T>(ry)) +
           v0 * (uu * y * math::invComplex(Complex<T>{T{1}, T{0}} + ry));
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

    const Complex<T> invbetaPlusD{math::invComplex(betaPlusD)};

    const Complex<T> betaMinusD{
        (std::real(beta * std::conj(D)) > T{0}) ? -sigma2uu * invbetaPlusD : beta - D
    };

    const Complex<T> DT{D * t};

    Complex<T> y;
    Complex<T> oneMinusEDT;

    if (std::norm(D) > std::numeric_limits<T>::epsilon() * (T{1} + std::abs(DT)))
    {
        y = math::expm1Complex(-DT) * T{0.5} * math::invComplex(D);
        oneMinusEDT = -T{2} * D * y;
    }
    else
    {
        y = tDivTwo;
        oneMinusEDT = DT;
    }

    const Complex<T> ry{-betaMinusD * y};

    const Complex<T> betaMinusDTimesT{betaMinusD * t};

    const Complex<T> A{
        kappaThetaDivSigma2 * (betaMinusDTimesT - T{2} * math::log1pComplex<T>(ry))
    };

    const Complex<T> B{uu * y * math::invComplex(T{1} + ry)};

    const Complex<T> g{betaMinusD * invbetaPlusD};

    const Complex<T> eDT = Complex<T>{T{1}, T{0}} - oneMinusEDT;

    const Complex<T> Q{T{1} - g * eDT};

    const Complex<T> R{T{1} - g};

    return CharFunCache<T>{
        .logPsi = A + v0 * B,
        .A = A,
        .B = B,
        .beta = beta,
        .D = D,
        .betaMinusD = betaMinusD,
        .uu = uu,
        .eDT = eDT,
        .oneMinusEDT = oneMinusEDT,
        .g = g,
        .Q = Q,
        .R = R,
        .S = betaMinusDTimesT - T{2} * math::log1pComplex((Q - R) * math::invComplex(R)),
        .denomG = betaPlusD * betaPlusD
    };
}

} // namespace uv::models::heston::price::detail

#include "Models/Heston/Price/Detail/CharFunction.inl"