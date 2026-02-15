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
#include "Math/Functions/Primitive.hpp"
#include "Models/Heston/Price/Detail/CharFunction.hpp"

#include <cmath>
#include <complex>

namespace uv::models::heston::price::detail
{
template <std::floating_point T> T Integrand<T>::operator()(T x) const noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};

    const Complex<T> h{iAlpha + x * onePlusITanPhi};

    const Complex<T> hMinusI{h - i};

    const Complex<T> logPsi{charFunction(
        kappa,
        kappaThetaDivSigma2,
        sigma2,
        v0,
        t,
        tDivTwo,
        sigmaRho,
        hMinusI
    )};

    const Complex<T> invDenom{math::invComplex(hMinusI * h)};

    return std::real(std::exp(logPsi + (x * c)) * onePlusITanPhi * invDenom);
}

template <std::floating_point T>
std::array<T, 6> BatchIntegrand<T>::operator()(T x) const noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};
    const Complex<T> h{iAlpha + x * onePlusITanPhi};
    const Complex<T> hMinusI{h - i};

    const auto cfData = charFunctionCached(
        kappa,
        kappaThetaDivSigma2,
        sigma2,
        v0,
        t,
        tDivTwo,
        sigmaRho,
        hMinusI
    );

    const Complex<T> betaMinusD{cfData.betaMinusD};
    const Complex<T> Q{cfData.Q};
    const Complex<T> invQ{math::invComplex(Q)};
    const Complex<T> eDT{cfData.eDT};
    const Complex<T> oneMinusEDT{cfData.oneMinusEDT};

    const GradientBD<T> gradientBD{
        .t = t,
        .invSigma2 = invSigma2,
        .invR = math::invComplex(cfData.R),
        .betaMinusDinvSigma2 = {betaMinusD * invSigma2},
        .deDTdD = {-t * eDT},
        .eDT = eDT,
        .oneMinusEDT = oneMinusEDT,
        .g = cfData.g,
        .invQ = invQ,
        .invQ2 = {invQ * invQ},
        .fracB = {oneMinusEDT * invQ},
        .Q = Q
    };

    const Complex<T> ui{hMinusI * i};

    const Complex<T> dbetaDs{-rho * ui};
    const Complex<T> dbetaDr{-sigma * ui};

    const Complex<T> d{cfData.D};
    const Complex<T> invD{math::invComplex(d)};
    const Complex<T> beta{cfData.beta};

    const Complex<T> dDdk{beta * invD};
    const Complex<T> dDds{(beta * dbetaDs + sigma * cfData.uu) * invD};
    const Complex<T> dDdr{dDdk * dbetaDr};

    const Complex<T> invDenomG{T{2} * math::invComplex(cfData.denomG)};

    const Complex<T> dgdk{(d * dbetaDk - beta * dDdk) * invDenomG};
    const Complex<T> dgds{(d * dbetaDs - beta * dDds) * invDenomG};
    const Complex<T> dgdr{(d * dbetaDr - beta * dDdr) * invDenomG};

    auto [dBdk, dSdk]{gradientBD.template eval<false>(dbetaDk, dDdk, dgdk)};
    auto [dBds, dSds]{
        gradientBD.template eval<true>(dbetaDs, dDds, dgds, betaMinusD * invSigma3Two)
    };
    auto [dBdr, dSdr]{gradientBD.template eval<false>(dbetaDr, dDdr, dgdr)};

    const Complex<T> s{cfData.S};

    const Complex<T> dAdk{dKdk * s + kappaThetaDivSigma2 * dSdk};
    const Complex<T> dAds{dKds * s + kappaThetaDivSigma2 * dSds};
    const Complex<T> dAdr{kappaThetaDivSigma2 * dSdr};

    const Complex<T> invDenom{math::invComplex(hMinusI * h)};

    const Complex<T> kernel{
        std::exp(cfData.logPsi + (x * c)) * onePlusITanPhi * invDenom
    };

    return {
        std::real(kernel),
        std::real(kernel * (dAdk + v0 * dBdk)),
        std::real(kernel * (cfData.A * invTheta)),
        std::real(kernel * (dAds + v0 * dBds)),
        std::real(kernel * (dAdr + v0 * dBdr)),
        std::real(kernel * cfData.B)
    };
}

template <std::floating_point T> template <bool HasSigmaTerm>
GradResult<T> GradientBD<T>::eval(
    const Complex<T> dbeta,
    const Complex<T> dD,
    const Complex<T> dg,
    const Complex<T> sigmaTerm
) const noexcept
{
    const Complex<T> a{-deDTdD * dD};
    const Complex<T> dbetaMinusdD{dbeta - dD};
    const Complex<T> commonGTerm{-dg * eDT + g * a};

    Complex<T> first{dbetaMinusdD * invSigma2};
    if constexpr (HasSigmaTerm)
        first += sigmaTerm;
    first *= fracB;

    return GradResult{
        first + betaMinusDinvSigma2 * ((a * Q - oneMinusEDT * commonGTerm) * invQ2),
        dbetaMinusdD * t - T{2} * (commonGTerm * invQ + dg * invR)
    };
}
} // namespace uv::models::heston::price::detail