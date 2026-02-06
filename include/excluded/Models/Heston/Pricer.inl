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

#include "Math/Functions.hpp"
#include "Utils/Aux/Errors.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace uv::models::heston
{
template <std::floating_point T, std::size_t N>
Pricer<T, N>::Pricer(
    std::shared_ptr<const math::integration::TanHSinH<T, N>> quad,
    const Config<T>& config
)
    : quad_(std::move(quad)),
      config_(config)
{

    UV_REQUIRE(
        (config_.alphaItm <= -T(1.0) - config_.eps) && (config_.alphaOtm >= config_.eps),
        ErrorCode::InvalidArgument,
        "Alpha must be outside [-1, 0] range"
    );
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::callPrice(T kappa, T theta, T sigma, T rho, T v0, T t, T F, T r, T K)
    const noexcept
{

    const T w{std::log(F / K)};

    const T alpha(getAlpha(w));

    const T R{Pricer<T, N>::getResidues(alpha, F, K)};

    const T phi{Pricer<T, N>::getPhi(kappa, theta, sigma, rho, v0, t, w)};

    constexpr Complex<T> i(T(0.0), T(1));
    const T tanPhi{std::tan(phi)};
    const Complex<T> onePlusITanPhi{T(1) + i * tanPhi};
    const Complex<T> c{(i - tanPhi) * w};
    const Complex<T> iAlpha{-i * alpha};

    auto integrand = [=](T x) noexcept -> T
    {
        const Complex<T> h{iAlpha + x * onePlusITanPhi};

        const Complex<T> hMinusI{h - i};

        const Complex<T> psi{
            Pricer<T, N>::charFunction(kappa, theta, sigma, rho, v0, t, hMinusI)
        };

        return std::real(std::exp(x * c) * psi / (hMinusI * h) * onePlusITanPhi);
    };

    constexpr T pi{std::numbers::pi_v<T>};
    return std::exp(-r * t) *
           (R - (F / pi) * std::exp(alpha * w) * quad_->integrateZeroToInf(integrand));
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::callPrice(T t, T F, T r, T K) const
{

    UV_REQUIRE(
        (params_.has_value()),
        ErrorCode::InvalidArgument,
        "Pricer::callPrice: Heston parameters not set. Call "
        "setParams(...) first."
    );

    const Params<T>& params{*params_};

    return callPrice(
        params.kappa,
        params.theta,
        params.sigma,
        params.rho,
        params.v0,
        t,
        F,
        r,
        K
    );
}

template <std::floating_point T, std::size_t N>
std::array<T, 6> Pricer<T, N>::callPriceWithGradient(
    T kappa,
    T theta,
    T sigma,
    T rho,
    T v0,
    T t,
    T F,
    T r,
    T K
) const noexcept
{

    const T w{std::log(F / K)};

    const T alpha(getAlpha(w));

    const T R{Pricer<T, N>::getResidues(alpha, F, K)};

    const T phi{Pricer<T, N>::getPhi(kappa, theta, sigma, rho, v0, t, w)};

    constexpr Complex<T> i(T(0), T(1));
    const T tanPhi{std::tan(phi)};
    const Complex<T> onePlusITanPhi{T(1) + i * tanPhi};
    const Complex<T> c{(i - tanPhi) * w};
    const Complex<T> iAlpha{-i * alpha};

    auto batchIntegrand = [=](T x) noexcept -> std::array<T, 6>
    {
        const Complex<T> h{iAlpha + x * onePlusITanPhi};
        const Complex<T> hMinusI{h - i};
        const Complex<T> kernel{std::exp(x * c) * onePlusITanPhi / (hMinusI * h)};

        const CharFunCache cfData{
            Pricer<T, N>::charFunctionCal(kappa, theta, sigma, rho, v0, t, hMinusI)
        };

        const Complex<T> psi{cfData.psi};
        const Complex<T> A{cfData.A};
        const Complex<T> B{cfData.B};
        const Complex<T> beta{cfData.beta};
        const Complex<T> D{cfData.D};
        const Complex<T> betaPlusD{cfData.betaPlusD};
        const Complex<T> betaMinusD{cfData.betaMinusD};
        const Complex<T> kFac{cfData.kFac};
        const Complex<T> ui{cfData.ui};
        const Complex<T> uu{cfData.uu};
        const Complex<T> eDT{cfData.eDT};
        const Complex<T> g{cfData.g};
        const Complex<T> Q{cfData.Q};
        const Complex<T> invQ{cfData.invQ};
        const Complex<T> invQ2{cfData.invQ2};
        const Complex<T> R{cfData.R};
        const Complex<T> S{cfData.S};
        const Complex<T> fracB{cfData.fracB};
        const Complex<T> denomG{cfData.denomG};
        const Complex<T> betaMinusDinvSigma2{cfData.betaMinusDinvSigma2};
        const T sigma2{cfData.sigma2};
        const T invSigma2{cfData.invSigma2};
        const T kappaTheta{cfData.kappaTheta};

        const T sigma3{sigma2 * sigma};
        const T invTheta{T(1) / theta};
        const Complex<T> u2{uu - ui};
        const Complex<T> deDT_dD{-t * eDT};

        constexpr Complex<T> dbeta_dk{T(1), T(0.0)};
        const Complex<T> dbeta_ds{-rho * ui};
        const Complex<T> dbeta_dr{-sigma * ui};

        const Complex<T> dD_dk{beta / D};
        const Complex<T> dD_ds{(beta * dbeta_ds + sigma * (ui + u2)) / D};
        const Complex<T> dD_dr{dD_dk * dbeta_dr};

        const auto dg_from = [&D, &denomG, &beta](
                                 const Complex<T>& dbeta,
                                 const Complex<T>& dD
                             ) noexcept -> Complex<T>
        {
            return T(2) * (D * dbeta - beta * dD) / denomG;
        };

        const Complex<T> dg_dk{dg_from(dbeta_dk, dD_dk)};
        const Complex<T> dg_ds{dg_from(dbeta_ds, dD_ds)};
        const Complex<T> dg_dr{dg_from(dbeta_dr, dD_dr)};

        const auto dB_from = [&](const Complex<T>& dbeta, const Complex<T>& dD, T dC
                             ) noexcept -> Complex<T>
        {
            const Complex<T> d_one_over_C2{(-T(2) * dC) / sigma3};
            const Complex<T> dpref{(dbeta - dD) * invSigma2 + betaMinusD * d_one_over_C2};

            const Complex<T> dN1{(+deDT_dD) * (-dD)};
            const Complex<T> dQ{-(dg_from(dbeta, dD) * eDT) - g * (deDT_dD * dD)};
            const Complex<T> dfrac{(dN1 * Q - (T(1) - eDT) * dQ) * invQ2};

            return dpref * fracB + betaMinusDinvSigma2 * dfrac;
        };

        constexpr T dC_dk{T(0.0)};
        constexpr T dC_ds{T(1)};
        constexpr T dC_dr{T(0.0)};

        const Complex<T> dB_dk{dB_from(dbeta_dk, dD_dk, dC_dk)};
        const Complex<T> dB_ds{dB_from(dbeta_ds, dD_ds, dC_ds)};
        const Complex<T> dB_dr{dB_from(dbeta_dr, dD_dr, dC_dr)};

        const Complex<T> dK_dk{theta * invSigma2};
        const Complex<T> dK_ds{(T(-T(2.0)) * kappaTheta) / sigma3};
        constexpr Complex<T> dK_dr{T(0.0), T(0.0)};

        const auto dS_from = [&](const Complex<T>& dbeta,
                                 const Complex<T>& dD,
                                 const Complex<T>& dg) noexcept -> Complex<T>
        {
            const Complex<T> dQ{-(dg * eDT + g * (deDT_dD * dD))};
            return (dbeta - dD) * t - T(2) * (dQ * invQ + dg / R);
        };

        const Complex<T> dS_dk{dS_from(dbeta_dk, dD_dk, dg_dk)};
        const Complex<T> dS_ds{dS_from(dbeta_ds, dD_ds, dg_ds)};
        const Complex<T> dS_dr{dS_from(dbeta_dr, dD_dr, dg_dr)};

        const Complex<T> dA_dk{dK_dk * S + kFac * dS_dk};
        const Complex<T> dA_ds{dK_ds * S + kFac * dS_ds};
        const Complex<T> dA_dr{dK_dr * S + kFac * dS_dr};

        const Complex<T> kernelPsi{kernel * psi};
        return {
            std::real(kernelPsi),
            std::real(kernelPsi * (dA_dk + v0 * dB_dk)),
            std::real(kernelPsi * (A * invTheta)),
            std::real(kernelPsi * (dA_ds + v0 * dB_ds)),
            std::real(kernelPsi * (dA_dr + v0 * dB_dr)),
            std::real(kernelPsi * B)
        };
    };

    const auto integrals = quad_->template integrateZeroToInfMulti<6>(batchIntegrand);

    const T disc{std::exp(-r * t)};
    const T pref{(F / std::numbers::pi_v<T>)*std::exp(alpha * w)};
    const T scale{disc * pref};

    std::array<T, 6> out{};
    out[0] = T(disc * (R - pref * integrals[0]));
    out[1] = T(-scale * integrals[1]);
    out[2] = T(-scale * integrals[2]);
    out[3] = T(-scale * integrals[3]);
    out[4] = T(-scale * integrals[4]);
    out[5] = T(-scale * integrals[5]);
    return out;
}

template <std::floating_point T, std::size_t N>
void Pricer<T, N>::setParams(const Params<T>& params) noexcept
{
    params_ = params;
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getResidues(T alpha, const T F, const T K) noexcept
{
    if (alpha < -T(1))
        return F - K;
    else
        return T(0.0);
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getAlpha(T w) const noexcept
{
    if (w >= T(0.0))
        return config_.alphaItm;
    else
        return config_.alphaOtm;
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::getPhi(T kappa, T theta, T sigma, T rho, T v0, T t, T w) noexcept
{
    if ((w * (rho - sigma * w / (v0 + kappa * theta * t))) >= T(0.0))
        return T(0.0);
    else
        return std::copysign(std::numbers::pi_v<T> / T(T(12.0)), w);
}

template <std::floating_point T, std::size_t N>
Complex<T> Pricer<T, N>::charFunction(
    T kappa,
    T theta,
    T sigma,
    T rho,
    T v0,
    T t,
    const Complex<T>& u
) noexcept
{

    constexpr Complex<T> i{T(0.0), T(1)};

    Complex<T> uu{u * (u + i)};

    const T sigma2{sigma * sigma};

    const Complex<T> beta{kappa - i * sigma * rho * u};

    const Complex<T> D{std::sqrt(beta * beta + sigma2 * uu)};

    const Complex<T> betaPlusD{beta + D};

    const Complex<T> r{
        (std::real(beta * std::conj(D)) > T(0.0)) ? -sigma2 * uu / betaPlusD : beta - D
    };

    const Complex<T> DT{D * t};
    Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T(1) + std::abs(DT)))
            ? math::expm1Complex(-DT) / (T(2) * D)
            : Complex<T>(-t / T(2))
    };

    const Complex<T> ry{r * y};

    const Complex<T> A{
        (kappa * theta / sigma2) * (r * t - T(2) * math::log1pComplex<T>(-ry))
    };

    const Complex<T> B{uu * y / (T(1) - ry)};

    return std::exp(A + v0 * B);
}

template <std::floating_point T, std::size_t N>
Pricer<T, N>::CharFunCache Pricer<T, N>::charFunctionCal(
    T kappa,
    T theta,
    T sigma,
    T rho,
    T v0,
    T t,
    const Complex<T>& u
) noexcept
{

    constexpr Complex<T> i{T(0.0), T(1)};

    const Complex<T> uu{u * (u + i)};

    const T sigma2{sigma * sigma};

    const T invSigma2{T(1) / sigma2};

    const Complex<T> ui{u * i};

    const Complex<T> beta{kappa - sigma * rho * ui};

    const Complex<T> D{std::sqrt(beta * beta + sigma2 * uu)};

    const Complex<T> betaPlusD{beta + D};

    const Complex<T> betaMinusD{
        (std::real(beta * std::conj(D)) > T(0.0)) ? -sigma2 * uu / betaPlusD : beta - D
    };

    const Complex<T> DT{D * t};

    const Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T(1) + std::abs(DT)))
            ? math::expm1Complex(-DT) / (T(2) * D)
            : Complex<T>(-t / T(2))
    };

    const Complex<T> ry{betaMinusD * y};

    const T kappaTheta{kappa * theta};

    const Complex<T> kFac{kappaTheta * invSigma2};

    const Complex<T> A{kFac * (betaMinusD * t - T(2) * math::log1pComplex<T>(-ry))};

    const Complex<T> B{uu * y / (T(1.0) - ry)};

    const Complex<T> eDT{std::exp(-DT)};

    const Complex<T> g{betaMinusD / betaPlusD};

    const Complex<T> Q{T(1) - g * eDT};

    const Complex<T> invQ{T(1) / Q};

    const Complex<T> R{T(1) - g};

    return CharFunCache{
        std::exp(A + v0 * B),
        A,
        B,
        beta,
        D,
        DT,
        betaPlusD,
        betaMinusD,
        ui,
        kFac,
        invSigma2,
        kappaTheta,
        sigma2,
        uu,
        eDT,
        betaMinusD / betaPlusD,
        Q,
        invQ,
        invQ * invQ,
        R,
        betaMinusD * t - T(2) * std::log(Q / R),
        (T(1) - eDT) * invQ,
        betaPlusD * betaPlusD,
        betaMinusD * invSigma2
    };
}
} // namespace uv::models::heston