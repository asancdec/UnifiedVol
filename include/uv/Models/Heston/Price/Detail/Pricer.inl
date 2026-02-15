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
    constexpr Complex<T> i{T{0}, T{1}};

    const T w{std::log(F / K)};
    const T alpha{getAlpha_(w)};
    const T tanPhi{std::tan(getPhi_(kappa, theta, sigma, rho, v0, t, w))};
    const T sigma2{sigma * sigma};
    const T invSigma3{T{1} / (sigma2 * sigma)};
    const T invSigma2{T{1} / sigma2};
    const T kappaTheta{kappa * theta};

    const detail::BatchIntegrand<T> batchIntegrand{
        .sigmaRho = {-i * sigma * rho},
        .tDivTwo = {-t * T{0.5}},
        .iAlpha = {-i * alpha},
        .onePlusITanPhi = {T{1} + i * tanPhi},
        .dbetaDk = {T{1}, T{0}},
        .c = {(i - tanPhi) * w},
        .kappa = kappa,
        .invSigma2 = invSigma2,
        .kappaThetaDivSigma2 = {kappaTheta / sigma2},
        .sigma = sigma,
        .sigma2 = sigma2,
        .rho = rho,
        .v0 = v0,
        .t = t,
        .invTheta = {T{1} / theta},
        .dKdk = {theta * invSigma2},
        .dKds = {T{-2} * kappaTheta * invSigma3},
        .invSigma3Two = {T{-2} * invSigma3}
    };

    const auto integrals = quad_->template integrateZeroToInfMulti<6>(batchIntegrand);

    constexpr T invPi{T{1} / std::numbers::pi_v<T>};
    const T pref{-(F * invPi) * std::exp(alpha * w)};
    const T scale{dF * pref};

    return std::array<T, 6>{
        dF * (getResidues_(alpha, F, K) + pref * integrals[0]),
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