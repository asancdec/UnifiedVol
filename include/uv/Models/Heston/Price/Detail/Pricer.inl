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

#include "Base/Macros/Require.hpp"
#include "Math/Functions/Primitive.hpp"
#include "Models/Heston/Price/Detail/Integrand.hpp"
#include <ceres/context.h>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace uv::models::heston::price
{

template <std::floating_point T, std::size_t N> Pricer<T, N>::Pricer()
    : Pricer(std::make_shared<const math::integration::TanHSinH<T, N>>())
{
    setAlphas_(Config<T>{});
    validateAlphaDomain_();
}

template <std::floating_point T, std::size_t N> Pricer<T, N>::Pricer(
    std::shared_ptr<const math::integration::TanHSinH<T, N>> quad,
    const Config<T>& config
)
    : quad_(std::move(quad))
{
    setAlphas_(config);
    UV_REQUIRE_NON_NULL(quad_);
    validateAlphaDomain_();
}

template <std::floating_point T, std::size_t N>
void Pricer<T, N>::Pricer::validateAlphaDomain_() const
{
    constexpr T EPS{std::numeric_limits<T>::epsilon() * 10};
    UV_REQUIRE_EQUAL_OR_LESS(alphaItm_, 1.0 - EPS);
    UV_REQUIRE_EQUAL_OR_GREATER(alphaOtm_, EPS);
}

template <std::floating_point T, std::size_t N>
void Pricer<T, N>::Pricer::validateCallPrice_(T t, T dF, T F, T K) const
{
    UV_REQUIRE_SET(params_);

    UV_REQUIRE_FINITE(t);
    UV_REQUIRE_FINITE(dF);
    UV_REQUIRE_FINITE(F);
    UV_REQUIRE_FINITE(K);

    UV_REQUIRE_POSITIVE(t);
    UV_REQUIRE_POSITIVE(dF);
}

template <std::floating_point T, std::size_t N>
void Pricer<T, N>::Pricer::setAlphas_(const Config<T>& config)
{
    alphaItm_ = config.alphaItm;
    alphaOtm_ = config.alphaOtm;
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getAlpha_(T w) const noexcept
{
    if (w >= T{0})
        return alphaItm_;

    return alphaOtm_;
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getResidues_(T alpha, const T F, const T K) noexcept
{
    if (alpha < -T{1})
        return F - K;

    return T{0};
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getPhi_(T kappa, T theta, T sigma, T rho, T v0, T t, T w) noexcept
{
    if (w * (rho - sigma * w / (v0 + kappa * theta * t)) >= T{0})
        return T{0};

    constexpr T piDivTwelve{std::numbers::pi_v<T> / T{12}};
    return std::copysign(piDivTwelve, w);
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::callPrice(T kappa, T theta, T sigma, T rho, T v0, T t, T dF, T F, T K)
    const noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};

    const T w{std::log(F / K)};
    const T alpha{getAlpha_(w)};
    const T tanPhi{std::tan(getPhi_(kappa, theta, sigma, rho, v0, t, w))};

    const T sigma2{sigma * sigma};

    const detail::Integrand<T> integrand{
        .iAlpha = {T{0}, -alpha},
        .onePlusITanPhi = {T{1}, tanPhi},
        .c = {-tanPhi * w, w},
        .tDivTwo = {-t * T{0.5}},
        .sigmaRho = {-i * (sigma * rho)},
        .kappa = kappa,
        .kappaThetaDivSigma2 = kappa * theta / sigma2,
        .sigma2 = sigma2,
        .v0 = v0,
        .t = t
    };

    constexpr T invPi{T{1} / std::numbers::pi_v<T>};

    return dF * (getResidues_(alpha, F, K) - (F * invPi) * std::exp(alpha * w) *
                                                 quad_->integrateZeroToInf(integrand));
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::callPrice(T t, T dF, T F, T K, bool doValidate) const
{
    if (doValidate)
        validateCallPrice_(t, dF, F, K);

    const Params<T>& params{*params_};

    return callPrice(
        params.kappa,
        params.theta,
        params.sigma,
        params.rho,
        params.v0,
        t,
        dF,
        F,
        K
    );
}

template <std::floating_point T, std::size_t N> void Pricer<T, N>::callPrice(
    std::span<T> out,
    T t,
    T dF,
    T F,
    std::span<const T> strikes,
    bool doValidate
) const
{
    if (doValidate)
    {
        UV_REQUIRE_NON_EMPTY(strikes);
        UV_REQUIRE_SAME_SIZE(out, strikes);
    }

    for (std::size_t i{0}; i < strikes.size(); ++i)
    {
        out[i] = callPrice(t, dF, F, strikes[i], doValidate);
    }
}

template <std::floating_point T, std::size_t N> core::Matrix<T> Pricer<T, N>::callPrice(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    bool doValidate
) const
{
    std::size_t numMaturities{volSurface.numMaturities()};

    std::span<const T> maturities{volSurface.maturities()};

    const Vector<T> discountFactors{curve.interpolateDF(maturities)};

    std::span<const T> forwards(volSurface.forwards());
    std::span<const T> strikes{volSurface.strikes()};

    core::Matrix<T> out{numMaturities, volSurface.numStrikes()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        callPrice(
            out[i],
            maturities[i],
            discountFactors[i],
            forwards[i],
            strikes,
            doValidate
        );
    }

    return out;
}

template <std::floating_point T, std::size_t N>
std::array<T, 6> Pricer<T, N>::callPriceWithGradient(
    T kappa,
    T theta,
    T sigma,
    T rho,
    T v0,
    T t,
    T dF,
    T F,
    T K
) const noexcept
{
    constexpr Complex<T> i(T{0}, T{1});

    const T w{std::log(F / K)};

    const T alpha(getAlpha_(w));
    const T R{getResidues_(alpha, F, K)};
    const T tanPhi{std::tan(getPhi_(kappa, theta, sigma, rho, v0, t, w))};
    const T sigma2{sigma * sigma};
    const T invSigma3{T{1} / (sigma2 * sigma)};
    const T invSigma2{T{1} / sigma2};
    const T invSigma3Two{T{-2} * invSigma3};
    const T kappaTheta{kappa * theta};
    const T invTheta{T{1} / theta};
    const T kappaThetaDivSigma2{kappaTheta / sigma2};
    const T dK_dk = theta * invSigma2;
    const T dK_ds = (T{-2} * kappaTheta) * invSigma3;

    const Complex<T> sigmaRho{-i * sigma * rho};
    const Complex<T> tDivTwo{-t * T{0.5}};
    const Complex<T> onePlusITanPhi{T{1} + i * tanPhi};
    const Complex<T> c{(i - tanPhi) * w};
    const Complex<T> iAlpha{-i * alpha};
    constexpr Complex<T> dbeta_dk{T{1}, T{0}};

    auto batchIntegrand = [kappa,
                           kappaTheta,
                           invSigma3,
                           sigmaRho,
                           invSigma2,
                           kappaThetaDivSigma2,
                           theta,
                           sigma,
                           sigma2,
                           rho,
                           v0,
                           t,
                           tDivTwo,
                           invTheta,
                           iAlpha,
                           onePlusITanPhi,
                           dK_dk,
                           dK_ds,
                           dbeta_dk,
                           invSigma3Two,
                           c](T x) noexcept -> std::array<T, 6>
    {
        constexpr Complex<T> i{T{0}, T{1}};

        const Complex<T> h{iAlpha + x * onePlusITanPhi};
        const Complex<T> hMinusI{h - i};
        const Complex<T> kernel{std::exp(x * c) * onePlusITanPhi / (hMinusI * h)};

        const detail::CharFunCache<T> cfData = detail::charFunctionCached(
            kappa,
            kappaThetaDivSigma2,
            sigma2,
            v0,
            t,
            tDivTwo,
            sigmaRho,
            hMinusI
        );

        const Complex<T> beta = cfData.beta;
        const Complex<T> D = cfData.D;
        const Complex<T> betaMinusD = cfData.betaMinusD;
        const Complex<T> eDT = cfData.eDT;
        const Complex<T> g = cfData.g;
        const Complex<T> invQ = cfData.invQ;
        const Complex<T> S = cfData.S;
        const Complex<T> fracB = cfData.fracB;
        const Complex<T> Q = cfData.Q;

        const Complex<T> ui = hMinusI * i;
        const Complex<T> u2 = cfData.uu - ui;
        const Complex<T> deDT_dD = -t * eDT;
        const Complex<T> oneMinusEDT{T{1} - eDT};
        const Complex<T> invD = T{1} / D;
        const Complex<T> invDenomG = T{2} / cfData.denomG;

        const Complex<T> dbeta_ds = -rho * ui;
        const Complex<T> dbeta_dr = -sigma * ui;

        const Complex<T> dD_dk = beta * invD;
        const Complex<T> dD_ds = (beta * dbeta_ds + sigma * (ui + u2)) * invD;
        const Complex<T> dD_dr = dD_dk * dbeta_dr;

        const Complex<T> dg_dk = (D * dbeta_dk - beta * dD_dk) * invDenomG;
        const Complex<T> dg_ds = (D * dbeta_ds - beta * dD_ds) * invDenomG;
        const Complex<T> dg_dr = (D * dbeta_dr - beta * dD_dr) * invDenomG;

        const Complex<T> betaMinusDinvSigma2 = betaMinusD * invSigma2;
        const Complex<T> invQ2 = invQ * invQ;

        const auto dB_from_0 = [invSigma2,
                                betaMinusDinvSigma2,
                                deDT_dD,
                                eDT,
                                oneMinusEDT,
                                g,
                                invQ2,
                                fracB,
                                Q](Complex<T> dbeta, Complex<T> dD, Complex<T> dg
                               ) noexcept -> Complex<T>
        {
            const Complex<T> a = -deDT_dD * dD;

            return (
                ((dbeta - dD) * invSigma2) * fracB +
                betaMinusDinvSigma2 *
                    ((a * Q - oneMinusEDT * (-(dg * eDT) + g * a)) * invQ2)
            );
        };
        const Complex<T> dB_dr = dB_from_0(dbeta_dr, dD_dr, dg_dr);
        const Complex<T> dB_dk = dB_from_0(dbeta_dk, dD_dk, dg_dk);

        // TODO: UNIFY dB and DS

        const Complex<T> a = -deDT_dD * dD_ds;
        const Complex<T> dB_ds =
            (((dbeta_ds - dD_ds) * invSigma2 + betaMinusD * invSigma3Two) * fracB +
             betaMinusDinvSigma2 *
                 ((a * Q - oneMinusEDT * (-(dg_ds * eDT) + g * a)) * invQ2));

        const auto dS_from = [eDT, g, deDT_dD, invQ, t, Rcf = cfData.R](
                                 Complex<T> dbeta,
                                 Complex<T> dD,
                                 Complex<T> dg
                             ) noexcept -> Complex<T>
        {
            const Complex<T> dQ = -(dg * eDT + g * (deDT_dD * dD));
            return (dbeta - dD) * t - T{2} * (dQ * invQ + dg / Rcf);
        };

        const Complex<T> dS_dk = dS_from(dbeta_dk, dD_dk, dg_dk);
        const Complex<T> dS_ds = dS_from(dbeta_ds, dD_ds, dg_ds);
        const Complex<T> dS_dr = dS_from(dbeta_dr, dD_dr, dg_dr);

        const Complex<T> dA_dk = dK_dk * S + kappaThetaDivSigma2 * dS_dk;
        const Complex<T> dA_ds = dK_ds * S + kappaThetaDivSigma2 * dS_ds;
        const Complex<T> dA_dr = kappaThetaDivSigma2 * dS_dr;

        const Complex<T> kernelPsi = kernel * cfData.psi;

        return {
            std::real(kernelPsi),
            std::real(kernelPsi * (dA_dk + v0 * dB_dk)),
            std::real(kernelPsi * (cfData.A * invTheta)),
            std::real(kernelPsi * (dA_ds + v0 * dB_ds)),
            std::real(kernelPsi * (dA_dr + v0 * dB_dr)),
            std::real(kernelPsi * cfData.B)
        };
    };

    const auto integrals = quad_->template integrateZeroToInfMulti<6>(batchIntegrand);

    constexpr T invPi{T{1} / std::numbers::pi_v<T>};
    const T pref{-(F * invPi) * std::exp(alpha * w)};
    const T scale{dF * pref};

    return std::array<T, 6>{
        dF * (R + pref * integrals[0]),
        scale * integrals[1],
        scale * integrals[2],
        scale * integrals[3],
        scale * integrals[4],
        scale * integrals[5]
    };
}

template <std::floating_point T, std::size_t N>
void Pricer<T, N>::setParams(const Params<T>& params) noexcept
{
    params_ = params;
}
} // namespace uv::models::heston::price