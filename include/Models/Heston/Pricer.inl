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


#include "Utils/Aux/Errors.hpp"    
#include "Math/Functions.hpp"

#include <utility>   
#include <cmath>     
#include <numbers>  
#include <limits>

namespace uv::models::heston
{
	template <std::size_t N>
	Pricer<N>::Pricer(std::shared_ptr<const math::TanHSinH<N>> quad,
		const Config& config) :
		quad_(std::move(quad)),
		config_(config)
	{
		// Alpha checks
		UV_REQUIRE
		(
			(config_.alphaItm <= -Real(1.0) - config_.eps) && (config_.alphaOtm >= config_.eps),
			ErrorCode::InvalidArgument,
			"Alpha must be outside [-1, 0] range"
		);
	}

	template <std::size_t N>
	Real Pricer<N>::callPrice(Real kappa,
		Real theta,
		Real sigma,
		Real rho,
		Real v0,
		Real T,
		Real F,
		Real r,
		Real K) const noexcept
	{
		// Calculate w
		const Real w{ std::log(F / K) };

		// Determine alpha
		const Real alpha(getAlpha(w));

		// Calculate residues
		const Real R{ Pricer<N>::getResidues(alpha, F, K) };

		// Determine phi
		const Real phi{ Pricer<N>::getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr Complex<Real> i(Real(Real(0.0)), Real(Real(1.0)));
		const Real tanPhi{ std::tan(phi) };
		const Complex<Real> onePlusITanPhi{ Real(Real(1.0)) + i * tanPhi };
		const Complex<Real> c{ (i - tanPhi) * w };
		const Complex<Real> iAlpha{- i * alpha};

		// Define the integrand
		auto integrand = [=](Real x) noexcept -> Real
			{
				// Calculate h(x)
				const Complex<Real> h{ iAlpha + x * onePlusITanPhi };

				// h - i
				const Complex<Real> hMinusI{ h - i };

				// Evaluate characteristic function at h(x) - i
				const Complex<Real> psi{ Pricer<N>::charFunction(kappa, theta, sigma, rho, v0, T, hMinusI) };

				// Calculate and return integrand
				return std::real(std::exp(x * c) * psi / (hMinusI * h) * onePlusITanPhi);
			};

		// Calculate and return call price
		constexpr Real pi{ std::numbers::pi_v<Real> };
		return std::exp(-r * T)
			* (R - (F / pi) * std::exp(alpha * w)
				* quad_->integrateZeroToInf(integrand));
	}

	template <std::size_t N>
	Real Pricer<N>::callPrice(Real T,
		Real F,
		Real r,
		Real K) const
	{
		// Will return error if class parameters are not set
		UV_REQUIRE
		(
			(params_.has_value()),
			ErrorCode::InvalidArgument,
			"Pricer::callPrice: Heston parameters not set. Call setParams(...) first."
		);

		// Dereference optional
		const Params& params{ *params_ };

		// Pass class instance parameters into the generic pricing function
		return callPrice(params.kappa,
			params.theta,
			params.sigma,
			params.rho,
			params.v0,
			T, F, r, K);
	}

	template <std::size_t N>
	std::array<Real, 6> Pricer<N>::callPriceWithGradient(Real kappa,
		Real theta,
		Real sigma,
		Real rho,
		Real v0,
		Real T,
		Real F,
		Real r,
		Real K) const noexcept
	{
		// Calculate w
		const Real w{ std::log(F / K) };

		// Determine alpha
		const Real alpha(getAlpha(w));

		// Calculate residues
		const Real R{ Pricer<N>::getResidues(alpha, F, K) };

		// Determine phi
		const Real phi{ Pricer<N>::getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr Complex<Real> i(Real(Real(0.0)), Real(Real(1.0)));
		const Real tanPhi{ std::tan(phi) };
		const Complex<Real> onePlusITanPhi{ Real(Real(1.0)) + i * tanPhi };
		const Complex<Real> c{ (i - tanPhi) * w };
		const Complex<Real> iAlpha{ -i * alpha };

		auto batchIntegrand = [=](Real x) noexcept -> std::array<Real, 6>
			{
				// Integrand precomputations
				const Complex<Real> h{ iAlpha + x * onePlusITanPhi };
				const Complex<Real> hMinusI{ h - i };
				const Complex<Real> kernel{ std::exp(x * c) * onePlusITanPhi / (hMinusI * h) };

				// Characteristic function and cached intermediates
				const CharFunData cfData
				{
					Pricer<N>::charFunctionCal(kappa, theta, sigma, rho, v0, T, hMinusI)
				};

				// Reuse intermediates
				const Complex<Real> psi{ cfData.psi };
				const Complex<Real> A{ cfData.A };
				const Complex<Real> B{ cfData.B };
				const Complex<Real> beta{ cfData.beta };
				const Complex<Real> D{ cfData.D };
				const Complex<Real> betaPlusD{ cfData.betaPlusD };
				const Complex<Real> betaMinusD{ cfData.betaMinusD };
				const Complex<Real> kFac{ cfData.kFac };
				const Complex<Real> ui{ cfData.ui };
				const Complex<Real> uu{ cfData.uu };
				const Complex<Real> eDT{ cfData.eDT };
				const Complex<Real> g{ cfData.g };
				const Complex<Real> Q{ cfData.Q };
				const Complex<Real> invQ{ cfData.invQ };
				const Complex<Real> invQ2{ cfData.invQ2 };
				const Complex<Real> R{ cfData.R };
				const Complex<Real> S{ cfData.S };
				const Complex<Real> fracB{ cfData.fracB };
				const Complex<Real> denomG{ cfData.denomG };
				const Complex<Real> betaMinusDinvSigma2{ cfData.betaMinusDinvSigma2 };
				const Real sigma2{ cfData.sigma2 };
				const Real invSigma2{ cfData.invSigma2 };
				const Real kappaTheta{ cfData.kappaTheta };

				// Precomputations
				const Real sigma3{ sigma2 * sigma };
				const Real invTheta{ Real(Real(1.0)) / theta };
				const Complex<Real> u2{ uu - ui };
				const Complex<Real> deDT_dD{ -T * eDT };

				// Derivatives of beta
				constexpr Complex<Real> dbeta_dk{ Real(Real(1.0)), Real(Real(0.0)) };
				const Complex<Real>     dbeta_ds{ -rho * ui };
				const Complex<Real>     dbeta_dr{ -sigma * ui };

				// Derivatives of D
				const Complex<Real> dD_dk{ beta / D };
				const Complex<Real> dD_ds{ (beta * dbeta_ds + sigma * (ui + u2)) / D };
				const Complex<Real> dD_dr{ dD_dk * dbeta_dr };

				// g' helper
				const auto dg_from = [&D, &denomG, &beta](const Complex<Real>& dbeta, const Complex<Real>& dD) noexcept -> Complex<Real>
					{
						return Real(Real(2.0)) * (D * dbeta - beta * dD) / denomG;
					};

				// Derivatives of g
				const Complex<Real> dg_dk{ dg_from(dbeta_dk, dD_dk) };
				const Complex<Real> dg_ds{ dg_from(dbeta_ds, dD_ds) };
				const Complex<Real> dg_dr{ dg_from(dbeta_dr, dD_dr) };

				// B' helper
				const auto dB_from = [&](const Complex<Real>& dbeta, const Complex<Real>& dD, Real dC) noexcept -> Complex<Real>
					{
						// Prefactor derivative wrt (beta, D, C=sigma)
						const Complex<Real> d_one_over_C2{ (-Real(Real(2.0)) * dC) / sigma3 };
						const Complex<Real> dpref{ (dbeta - dD) * invSigma2 + betaMinusD * d_one_over_C2 };

						// Fraction derivative using cached inverses
						const Complex<Real> dN1{ (+deDT_dD) * (-dD) };
						const Complex<Real> dQ{ -(dg_from(dbeta, dD) * eDT) - g * (deDT_dD * dD) };
						const Complex<Real> dfrac{ (dN1 * Q - (Real(Real(1.0)) - eDT) * dQ) * invQ2 };

						return dpref * fracB + betaMinusDinvSigma2 * dfrac;
					};

				// dC for each parameter (C := sigma)
				constexpr Real dC_dk{ Real(Real(0.0)) };
				constexpr Real dC_ds{ Real(Real(1.0)) };
				constexpr Real dC_dr{ Real(Real(0.0)) };

				// B partials
				const Complex<Real> dB_dk{ dB_from(dbeta_dk, dD_dk, dC_dk) };
				const Complex<Real> dB_ds{ dB_from(dbeta_ds, dD_ds, dC_ds) };
				const Complex<Real> dB_dr{ dB_from(dbeta_dr, dD_dr, dC_dr) };

				// A' pieces
				const Complex<Real> dK_dk{ theta * invSigma2 };
				const Complex<Real> dK_ds{ (Real(-Real(2.0)) * kappaTheta) / sigma3 };
				constexpr Complex<Real> dK_dr{ Real(Real(0.0)), Real(Real(0.0))};

				// S' helper
				const auto dS_from = [&](const Complex<Real>& dbeta, const Complex<Real>& dD, const Complex<Real>& dg) noexcept -> Complex<Real>
					{
						const Complex<Real> dQ{ -(dg * eDT + g * (deDT_dD * dD)) };
						return (dbeta - dD) * T - Real(Real(2.0)) * (dQ * invQ + dg / R);
					};

				// S partials
				const Complex<Real> dS_dk{ dS_from(dbeta_dk, dD_dk, dg_dk) };
				const Complex<Real> dS_ds{ dS_from(dbeta_ds, dD_ds, dg_ds) };
				const Complex<Real> dS_dr{ dS_from(dbeta_dr, dD_dr, dg_dr) };

				// A partials
				const Complex<Real> dA_dk{ dK_dk * S + kFac * dS_dk };
				const Complex<Real> dA_ds{ dK_ds * S + kFac * dS_ds };
				const Complex<Real> dA_dr{ dK_dr * S + kFac * dS_dr };

				// Outputs (real parts)
				const Complex<Real> kernelPsi{ kernel * psi };
				return
				{
					std::real(kernelPsi),                          // price integrand
					std::real(kernelPsi * (dA_dk + v0 * dB_dk)),   // dP/dkappa
					std::real(kernelPsi * (A * invTheta)),         // dP/dtheta
					std::real(kernelPsi * (dA_ds + v0 * dB_ds)),   // dP/dsigma
					std::real(kernelPsi * (dA_dr + v0 * dB_dr)),   // dP/drho
					std::real(kernelPsi * B)                       // dP/dv0
				};
			};

		// Integrate all 6 components in one loop
		const auto integrals = quad_->template integrateZeroToInfMulti<6>(batchIntegrand);

		// Precomputations
		const Real disc{ std::exp(-r * T) };
		const Real pref{ (F / std::numbers::pi_v<Real>) * std::exp(alpha * w) };
		const Real scale{ disc * pref };

		// Assemble price and gradients (unrolled)
		std::array<Real, 6> out{};
		out[0] = Real(disc * (R - pref * integrals[0]));
		out[1] = Real(-scale * integrals[1]);
		out[2] = Real(-scale * integrals[2]);
		out[3] = Real(-scale * integrals[3]);
		out[4] = Real(-scale * integrals[4]);
		out[5] = Real(-scale * integrals[5]);
		return out;
	}

	template <std::size_t N>
	template <typename T>
	void Pricer<N>::setParams(const std::array<T, 5>& params) noexcept
	{
		params_ = Params
		{
			Real(params[0]), // kappa
			Real(params[1]), // theta
			Real(params[2]), // sigma
			Real(params[3]), // rho
			Real(params[4])  // v0 
		};
	}

	template <std::size_t N>
	Real Pricer<N>::getResidues(Real alpha,
		const Real F,
		const Real K) noexcept
	{
		if (alpha < -Real(Real(1.0))) return F - K;
		else return Real(Real(0.0));
	}

	template <std::size_t N>
	Real Pricer<N>::getAlpha(Real w) const noexcept
	{
		if (w >= Real(Real(0.0))) return config_.alphaItm;
		else return config_.alphaOtm;
	}

	template <std::size_t N>
	Real Pricer<N>::getPhi(Real kappa,
		Real theta,
		Real sigma,
		Real rho,
		Real v0,
		Real T,
		Real w) noexcept
	{
		if ((w * (rho - sigma * w / (v0 + kappa * theta * T))) >= Real(Real(0.0))) return Real(Real(0.0));
		else return std::copysign(std::numbers::pi_v<Real> / Real(Real(12.0)), w);
	}

	template <std::size_t N>
	Complex<Real> Pricer<N>::charFunction(Real kappa,
		Real theta,
		Real sigma,
		Real rho,
		Real v0,
		Real T,
		const Complex<Real>& u) noexcept
	{
		// Define i 
		constexpr Complex<Real> i{ Real(Real(0.0)), Real(Real(1.0)) };

		// u * ( u+ 1)
		Complex<Real> uu{ u * (u + i) };

		// sigma^2
		const Real sigma2{ sigma * sigma };

		// beta := kappa - i * sigma * rho * u
		const Complex<Real> beta{ kappa - i * sigma * rho * u };

		// D : = sqrt(beta^2 + sigma^2 * u * (u + i))
		const Complex<Real> D{ std::sqrt(beta * beta + sigma2 * uu) };

		// beta + D
		const Complex<Real> betaPlusD{ beta + D };

		// beta - D
		const Complex<Real> r
		{
			(std::real(beta * std::conj(D)) > Real(Real(0.0)))
			? -sigma2 * uu / betaPlusD
			: beta - D
		};

		// y := expm1(−D * T) / (2 * D) with D≈0 fallback
		const Complex<Real> DT{ D * T };
		Complex<Real> y
		{
			(std::norm(D) > std::numeric_limits<Real>::epsilon() * (Real(Real(1.0)) + std::abs(DT)))
			? math::expm1Complex(-DT) / (Real(Real(2.0)) * D)
			: Complex<Real>(-T / Real(Real(2.0)))
		};

		// r * y 
		const Complex<Real> ry{ r * y };

		// A := (κ * theta / sigma^2) * (r * T − 2 * log1p(−r * y))
		const Complex<Real> A{(kappa * theta / sigma2) * (r * T - Real(Real(2.0)) * math::log1pComplex<Real>(-ry))	};

		// B := u * (u + i) * y / (1 − r * y)
		const Complex<Real> B{ uu * y / (Real(Real(1.0)) - ry) };

		// psi(u) := e^( A + v0 * B),
		return std::exp(A + v0 * B);
	}

	template <std::size_t N>
	CharFunData Pricer<N>::charFunctionCal(Real kappa,
		Real theta,
		Real sigma,
		Real rho,
		Real v0,
		Real T,
		const Complex<Real>& u) noexcept
	{
		// i
		constexpr Complex<Real> i{ Real(Real(0.0)), Real(Real(1.0)) };

		// u * (u + i)
		const Complex<Real> uu{ u * (u + i) };

		// sigma^2
		const Real sigma2{ sigma * sigma };

		// 1 / sigma^2
		const Real invSigma2{ Real(Real(1.0)) / sigma2 };

		// u * i
		const Complex<Real> ui{ u * i };

		// beta := kappa - sigma * rho * i*u
		const Complex<Real> beta{ kappa - sigma * rho * ui };

		// D := sqrt(beta^2 + sigma^2 * u(u+i))
		const Complex<Real> D{ std::sqrt(beta * beta + sigma2 * uu) };

		// beta + D
		const Complex<Real> betaPlusD{ beta + D };

		// beta - D  (stable branch)
		const Complex<Real> betaMinusD
		{
			(std::real(beta * std::conj(D)) > Real(Real(0.0)))
			? -sigma2 * uu / betaPlusD
			: beta - D
		};

		// DT := D * T
		const Complex<Real> DT{ D * T };

		// y := expm1(-DT) / (2D)   with D≈0 fallback
		const Complex<Real> y
		{
			(std::norm(D) > std::numeric_limits<Real>::epsilon() *
							  (Real(Real(1.0)) + std::abs(DT)))
			? math::expm1Complex(-DT) / (Real(Real(2.0)) * D)
			: Complex<Real>(-T / Real(Real(2.0)))
		};

		// r * y
		const Complex<Real> ry{ betaMinusD * y };

		// kappa * theta
		const Real kappaTheta{ kappa * theta };

		// (kappa * theta) / sigma^2
		const Complex<Real> kFac{ kappaTheta * invSigma2 };

		// A := (kappa*theta/sigma^2) * (r*T − 2 log(1 - r*y))
		const Complex<Real> A{ kFac * (betaMinusD * T - Real(Real(2.0))* math::log1pComplex<Real>(-ry)) };

		// B := u(u+i)*y / (1 − r*y)
		const Complex<Real> B{ uu * y / (Real(Real(1.0)) - ry) };

		// Rescued intermediates for gradient path

		// exp(-DT)
		const Complex<Real> eDT{ std::exp(-DT) };

		// g := (betaMinusD)/(betaPlusD)
		const Complex<Real> g{ betaMinusD / betaPlusD };

		// Q := 1 - g * eDT
		const Complex<Real> Q{ Real(Real(1.0)) - g * eDT };

		// 1/Q and 1/Q^2
		const Complex<Real> invQ{ Real(Real(1.0)) / Q };

		// R := 1 - g
		const Complex<Real> R{ Real(Real(1.0)) - g };

		// Return full CharFunData
		return
		{
			std::exp(A + v0 * B),									 // psi(u) := exp( A + v0 * B )
			A,													 	 // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
			B,														 // B := u(u+i)*y / (1 − r*y)
			beta,													 // beta := kappa − sigma*rho*(i*u)
			D,														 // D := sqrt( beta^2 + sigma^2*u(u+i) )
			DT,														 // DT := D*T
			betaPlusD,												// beta + D
			betaMinusD,												// beta − D   (stable)
			ui,														// ui := u*i
			kFac,													// kFac := (kappa*theta)/sigma^2
			invSigma2,												// 1 / sigma^2
			kappaTheta,												// kappa * theta
			sigma2,                                     // sigma^2
			uu,                                         // uu := u(u+i)
			eDT,                                        // eDT := exp( −DT )
			betaMinusD / betaPlusD,                     // g := (beta−D)/(beta+D)
			Q,                                          // Q := 1 − g*eDT
			invQ,                                       // 1 / Q
			invQ * invQ,                                // 1 / Q^2
			R,                                          // R := 1 − g
			betaMinusD * T - Real(Real(2.0)) * std::log(Q / R),    // S := (beta−D)*T − 2*log(Q/R)
			(Real(Real(1.0)) - eDT) * invQ,                        // fracB := (1 − eDT) / Q
			betaPlusD * betaPlusD,                      // denomG := (beta + D)^2
			betaMinusD * invSigma2                      // betaMinusD / sigma^2
		};
	}
}