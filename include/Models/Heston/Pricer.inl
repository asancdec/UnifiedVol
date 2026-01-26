// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.inl
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
    // Alpha checks
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
    // Calculate w
    const T w{std::log(F / K)};

    // Determine alpha
    const T alpha(getAlpha(w));

    // Calculate residues
    const T R{Pricer<T, N>::getResidues(alpha, F, K)};

    // Determine phi
    const T phi{Pricer<T, N>::getPhi(kappa, theta, sigma, rho, v0, t, w)};

    // Precomputations
    constexpr Complex<T> i(T(0.0), T(1));
    const T tanPhi{std::tan(phi)};
    const Complex<T> onePlusITanPhi{T(1) + i * tanPhi};
    const Complex<T> c{(i - tanPhi) * w};
    const Complex<T> iAlpha{-i * alpha};

    // Define the integrand
    auto integrand = [=](T x) noexcept -> T
    {
        // Calculate h(x)
        const Complex<T> h{iAlpha + x * onePlusITanPhi};

        // h - i
        const Complex<T> hMinusI{h - i};

        // Evaluate characteristic function at h(x) - i
        const Complex<T> psi{
            Pricer<T, N>::charFunction(kappa, theta, sigma, rho, v0, t, hMinusI)
        };

        // Calculate and return integrand
        return std::real(std::exp(x * c) * psi / (hMinusI * h) * onePlusITanPhi);
    };

    // Calculate and return call price
    constexpr T pi{std::numbers::pi_v<T>};
    return std::exp(-r * t) *
           (R - (F / pi) * std::exp(alpha * w) * quad_->integrateZeroToInf(integrand));
}

template <std::floating_point T, std::size_t N>
T Pricer<T, N>::callPrice(T t, T F, T r, T K) const
{
    // Will return error if class parameters are not set
    UV_REQUIRE(
        (params_.has_value()),
        ErrorCode::InvalidArgument,
        "Pricer::callPrice: Heston parameters not set. Call "
        "setParams(...) first."
    );

    // Dereference optional
    const Params<T>& params{*params_};

    // Pass class instance parameters into the generic pricing function
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
    // Calculate w
    const T w{std::log(F / K)};

    // Determine alpha
    const T alpha(getAlpha(w));

    // Calculate residues
    const T R{Pricer<T, N>::getResidues(alpha, F, K)};

    // Determine phi
    const T phi{Pricer<T, N>::getPhi(kappa, theta, sigma, rho, v0, t, w)};

    // Precomputations
    constexpr Complex<T> i(T(0), T(1));
    const T tanPhi{std::tan(phi)};
    const Complex<T> onePlusITanPhi{T(1) + i * tanPhi};
    const Complex<T> c{(i - tanPhi) * w};
    const Complex<T> iAlpha{-i * alpha};

    auto batchIntegrand = [=](T x) noexcept -> std::array<T, 6>
    {
        // Integrand precomputations
        const Complex<T> h{iAlpha + x * onePlusITanPhi};
        const Complex<T> hMinusI{h - i};
        const Complex<T> kernel{std::exp(x * c) * onePlusITanPhi / (hMinusI * h)};

        // Characteristic function and cached intermediates
        const CharFunCache cfData{
            Pricer<T, N>::charFunctionCal(kappa, theta, sigma, rho, v0, t, hMinusI)
        };

        // Reuse intermediates
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

        // Precomputations
        const T sigma3{sigma2 * sigma};
        const T invTheta{T(1) / theta};
        const Complex<T> u2{uu - ui};
        const Complex<T> deDT_dD{-t * eDT};

        // Derivatives of beta
        constexpr Complex<T> dbeta_dk{T(1), T(0.0)};
        const Complex<T> dbeta_ds{-rho * ui};
        const Complex<T> dbeta_dr{-sigma * ui};

        // Derivatives of D
        const Complex<T> dD_dk{beta / D};
        const Complex<T> dD_ds{(beta * dbeta_ds + sigma * (ui + u2)) / D};
        const Complex<T> dD_dr{dD_dk * dbeta_dr};

        // g' helper
        const auto dg_from = [&D, &denomG, &beta](
                                 const Complex<T>& dbeta,
                                 const Complex<T>& dD
                             ) noexcept -> Complex<T>
        {
            return T(2) * (D * dbeta - beta * dD) / denomG;
        };

        // Derivatives of g
        const Complex<T> dg_dk{dg_from(dbeta_dk, dD_dk)};
        const Complex<T> dg_ds{dg_from(dbeta_ds, dD_ds)};
        const Complex<T> dg_dr{dg_from(dbeta_dr, dD_dr)};

        // B' helper
        const auto dB_from = [&](const Complex<T>& dbeta, const Complex<T>& dD, T dC
                             ) noexcept -> Complex<T>
        {
            // Prefactor derivative wrt (beta, D, C=sigma)
            const Complex<T> d_one_over_C2{(-T(2) * dC) / sigma3};
            const Complex<T> dpref{(dbeta - dD) * invSigma2 + betaMinusD * d_one_over_C2};

            // Fraction derivative using cached inverses
            const Complex<T> dN1{(+deDT_dD) * (-dD)};
            const Complex<T> dQ{-(dg_from(dbeta, dD) * eDT) - g * (deDT_dD * dD)};
            const Complex<T> dfrac{(dN1 * Q - (T(1) - eDT) * dQ) * invQ2};

            return dpref * fracB + betaMinusDinvSigma2 * dfrac;
        };

        // dC for each parameter (C := sigma)
        constexpr T dC_dk{T(0.0)};
        constexpr T dC_ds{T(1)};
        constexpr T dC_dr{T(0.0)};

        // B partials
        const Complex<T> dB_dk{dB_from(dbeta_dk, dD_dk, dC_dk)};
        const Complex<T> dB_ds{dB_from(dbeta_ds, dD_ds, dC_ds)};
        const Complex<T> dB_dr{dB_from(dbeta_dr, dD_dr, dC_dr)};

        // A' pieces
        const Complex<T> dK_dk{theta * invSigma2};
        const Complex<T> dK_ds{(T(-T(2.0)) * kappaTheta) / sigma3};
        constexpr Complex<T> dK_dr{T(0.0), T(0.0)};

        // S' helper
        const auto dS_from = [&](const Complex<T>& dbeta,
                                 const Complex<T>& dD,
                                 const Complex<T>& dg) noexcept -> Complex<T>
        {
            const Complex<T> dQ{-(dg * eDT + g * (deDT_dD * dD))};
            return (dbeta - dD) * t - T(2) * (dQ * invQ + dg / R);
        };

        // S partials
        const Complex<T> dS_dk{dS_from(dbeta_dk, dD_dk, dg_dk)};
        const Complex<T> dS_ds{dS_from(dbeta_ds, dD_ds, dg_ds)};
        const Complex<T> dS_dr{dS_from(dbeta_dr, dD_dr, dg_dr)};

        // A partials
        const Complex<T> dA_dk{dK_dk * S + kFac * dS_dk};
        const Complex<T> dA_ds{dK_ds * S + kFac * dS_ds};
        const Complex<T> dA_dr{dK_dr * S + kFac * dS_dr};

        // Outputs (real parts)
        const Complex<T> kernelPsi{kernel * psi};
        return {
            std::real(kernelPsi),                        // price integrand
            std::real(kernelPsi * (dA_dk + v0 * dB_dk)), // dP/dkappa
            std::real(kernelPsi * (A * invTheta)),       // dP/dtheta
            std::real(kernelPsi * (dA_ds + v0 * dB_ds)), // dP/dsigma
            std::real(kernelPsi * (dA_dr + v0 * dB_dr)), // dP/drho
            std::real(kernelPsi * B)                     // dP/dv0
        };
    };

    // Integrate all 6 components in one loop
    const auto integrals = quad_->template integrateZeroToInfMulti<6>(batchIntegrand);

    // Precomputations
    const T disc{std::exp(-r * t)};
    const T pref{(F / std::numbers::pi_v<T>)*std::exp(alpha * w)};
    const T scale{disc * pref};

    // Assemble price and gradients (unrolled)
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
    // Define i
    constexpr Complex<T> i{T(0.0), T(1)};

    // u * ( u+ 1)
    Complex<T> uu{u * (u + i)};

    // sigma^2
    const T sigma2{sigma * sigma};

    // beta := kappa - i * sigma * rho * u
    const Complex<T> beta{kappa - i * sigma * rho * u};

    // D : = sqrt(beta^2 + sigma^2 * u * (u + i))
    const Complex<T> D{std::sqrt(beta * beta + sigma2 * uu)};

    // beta + D
    const Complex<T> betaPlusD{beta + D};

    // beta - D
    const Complex<T> r{
        (std::real(beta * std::conj(D)) > T(0.0)) ? -sigma2 * uu / betaPlusD : beta - D
    };

    // y := expm1(−D * T) / (2 * D) with D≈0 fallback
    const Complex<T> DT{D * t};
    Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T(1) + std::abs(DT)))
            ? math::expm1Complex(-DT) / (T(2) * D)
            : Complex<T>(-t / T(2))
    };

    // r * y
    const Complex<T> ry{r * y};

    // A := (κ * theta / sigma^2) * (r * T − 2 * log1p(−r * y))
    const Complex<T> A{
        (kappa * theta / sigma2) * (r * t - T(2) * math::log1pComplex<T>(-ry))
    };

    // B := u * (u + i) * y / (1 − r * y)
    const Complex<T> B{uu * y / (T(1) - ry)};

    // psi(u) := e^( A + v0 * B),
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
    // i
    constexpr Complex<T> i{T(0.0), T(1)};

    // u * (u + i)
    const Complex<T> uu{u * (u + i)};

    // sigma^2
    const T sigma2{sigma * sigma};

    // 1 / sigma^2
    const T invSigma2{T(1) / sigma2};

    // u * i
    const Complex<T> ui{u * i};

    // beta := kappa - sigma * rho * i*u
    const Complex<T> beta{kappa - sigma * rho * ui};

    // D := sqrt(beta^2 + sigma^2 * u(u+i))
    const Complex<T> D{std::sqrt(beta * beta + sigma2 * uu)};

    // beta + D
    const Complex<T> betaPlusD{beta + D};

    // beta - D  (stable branch)
    const Complex<T> betaMinusD{
        (std::real(beta * std::conj(D)) > T(0.0)) ? -sigma2 * uu / betaPlusD : beta - D
    };

    // DT := D * T
    const Complex<T> DT{D * t};

    // y := expm1(-DT) / (2D)   with D≈0 fallback
    const Complex<T> y{
        (std::norm(D) > std::numeric_limits<T>::epsilon() * (T(1) + std::abs(DT)))
            ? math::expm1Complex(-DT) / (T(2) * D)
            : Complex<T>(-t / T(2))
    };

    // r * y
    const Complex<T> ry{betaMinusD * y};

    // kappa * theta
    const T kappaTheta{kappa * theta};

    // (kappa * theta) / sigma^2
    const Complex<T> kFac{kappaTheta * invSigma2};

    // A := (kappa*theta/sigma^2) * (r*T − 2 log(1 - r*y))
    const Complex<T> A{kFac * (betaMinusD * t - T(2) * math::log1pComplex<T>(-ry))};

    // B := u(u+i)*y / (1 − r*y)
    const Complex<T> B{uu * y / (T(1.0) - ry)};

    // Rescued intermediates for gradient path

    // exp(-DT)
    const Complex<T> eDT{std::exp(-DT)};

    // g := (betaMinusD)/(betaPlusD)
    const Complex<T> g{betaMinusD / betaPlusD};

    // Q := 1 - g * eDT
    const Complex<T> Q{T(1) - g * eDT};

    // 1/Q and 1/Q^2
    const Complex<T> invQ{T(1) / Q};

    // R := 1 - g
    const Complex<T> R{T(1) - g};

    return CharFunCache{
        std::exp(A + v0 * B),   // psi(u) := exp( A + v0 * B )
        A,                      // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
        B,                      // B := u(u+i)*y / (1 − r*y)
        beta,                   // beta := kappa − sigma*rho*(i*u)
        D,                      // D := sqrt( beta^2 + sigma^2*u(u+i) )
        DT,                     // DT := D*T
        betaPlusD,              // beta + D
        betaMinusD,             // beta − D   (stable)
        ui,                     // ui := u*i
        kFac,                   // kFac := (kappa*theta)/sigma^2
        invSigma2,              // 1 / sigma^2
        kappaTheta,             // kappa * theta
        sigma2,                 // sigma^2
        uu,                     // uu := u(u+i)
        eDT,                    // eDT := exp( −DT )
        betaMinusD / betaPlusD, // g := (beta−D)/(beta+D)
        Q,                      // Q := 1 − g*eDT
        invQ,                   // 1 / Q
        invQ * invQ,            // 1 / Q^2
        R,                      // R := 1 − g
        betaMinusD * t - T(2) * std::log(Q / R), // S := (beta−D)*T − 2*log(Q/R)
        (T(1) - eDT) * invQ,                     // fracB := (1 − eDT) / Q
        betaPlusD * betaPlusD,                   // denomG := (beta + D)^2
        betaMinusD * invSigma2                   // betaMinusD / sigma^2
    };
}
} // namespace uv::models::heston