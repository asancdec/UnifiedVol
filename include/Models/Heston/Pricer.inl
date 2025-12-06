/**
* Pricer.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Aux/Errors.hpp"    
#include "Math/Functions.hpp"

#include <utility>   
#include <complex>  
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
			(config_.alphaItm <= -1.0 - config_.eps) && (config_.alphaOtm >= config_.eps),
			ErrorCode::InvalidArgument,
			"Alpha must be outside [-1, 0] range"
		);
	}

	template <std::size_t N>
	double Pricer<N>::callPrice(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double F,
		long double r,
		long double K) const noexcept
	{
		// Calculate w
		const long double w{ std::log(F / K) };

		// Determine alpha
		const long double alpha(getAlpha(w));

		// Calculate residues
		const long double R{ Pricer<N>::getResidues(alpha, F, K) };

		// Determine phi
		const long double phi{ Pricer<N>::getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr cplx i(0.0L, 1.0L);
		const long double tanPhi{ std::tan(phi) };
		const cplx onePlusITanPhi{ 1.0L + i * tanPhi };
		const cplx c{ (i - tanPhi) * w };
		const cplx iAlpha{- i * alpha};

		// Define the integrand
		auto integrand = [=](long double x) noexcept -> long double
			{
				// Calculate h(x)
				const cplx h{ iAlpha + x * onePlusITanPhi };

				// h - i
				const cplx hMinusI{ h - i };

				// Evaluate characteristic function at h(x) - i
				const cplx psi{ Pricer<N>::charFunction(kappa, theta, sigma, rho, v0, T, hMinusI) };

				// Calculate and return integrand
				return std::real(std::exp(x * c) * psi / (hMinusI * h) * onePlusITanPhi);
			};

		// Calculate and return call price
		constexpr long double pi{ std::numbers::pi_v<long double> };
		return std::exp(-r * T)
			* (R - (F / pi) * std::exp(alpha * w)
				* quad_->integrateZeroToInf(integrand));
	}

	template <std::size_t N>
	double Pricer<N>::callPrice(long double T,
		long double F,
		long double r,
		long double K) const
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
	std::array<double, 6> Pricer<N>::callPriceWithGradient(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double F,
		long double r,
		long double K) const noexcept
	{
		// Calculate w
		const long double w{ std::log(F / K) };

		// Determine alpha
		const long double alpha(getAlpha(w));

		// Calculate residues
		const long double R{ Pricer<N>::getResidues(alpha, F, K) };

		// Determine phi
		const long double phi{ Pricer<N>::getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr cplx i(0.0L, 1.0L);
		const long double tanPhi{ std::tan(phi) };
		const cplx onePlusITanPhi{ 1.0L + i * tanPhi };
		const cplx c{ (i - tanPhi) * w };
		const cplx iAlpha{ -i * alpha };

		auto batchIntegrand = [=](long double x) noexcept -> std::array<long double, 6>
			{
				// Integrand precomputations
				const cplx h{ iAlpha + x * onePlusITanPhi };
				const cplx hMinusI{ h - i };
				const cplx kernel{ std::exp(x * c) * onePlusITanPhi / (hMinusI * h) };

				// Characteristic function and cached intermediates
				const CharFunData cfData
				{
					Pricer<N>::charFunctionCal(kappa, theta, sigma, rho, v0, T, hMinusI)
				};

				// Reuse intermediates
				const cplx psi{ cfData.psi };
				const cplx A{ cfData.A };
				const cplx B{ cfData.B };
				const cplx beta{ cfData.beta };
				const cplx D{ cfData.D };
				const cplx betaPlusD{ cfData.betaPlusD };
				const cplx betaMinusD{ cfData.betaMinusD };
				const cplx kFac{ cfData.kFac };
				const cplx ui{ cfData.ui };
				const cplx uu{ cfData.uu };
				const cplx eDT{ cfData.eDT };
				const cplx g{ cfData.g };
				const cplx Q{ cfData.Q };
				const cplx invQ{ cfData.invQ };
				const cplx invQ2{ cfData.invQ2 };
				const cplx R{ cfData.R };
				const cplx S{ cfData.S };
				const cplx fracB{ cfData.fracB };
				const cplx denomG{ cfData.denomG };
				const cplx betaMinusDinvSigma2{ cfData.betaMinusDinvSigma2 };
				const long double sigma2{ cfData.sigma2 };
				const long double invSigma2{ cfData.invSigma2 };
				const long double kappaTheta{ cfData.kappaTheta };

				// Precomputations
				const long double sigma3{ sigma2 * sigma };
				const long double invTheta{ 1.0L / theta };
				const cplx u2{ uu - ui };
				const cplx deDT_dD{ -T * eDT };

				// Derivatives of beta
				constexpr cplx dbeta_dk{ 1.0L, 0.0L };
				const cplx     dbeta_ds{ -rho * ui };
				const cplx     dbeta_dr{ -sigma * ui };

				// Derivatives of D
				const cplx dD_dk{ beta / D };
				const cplx dD_ds{ (beta * dbeta_ds + sigma * (ui + u2)) / D };
				const cplx dD_dr{ dD_dk * dbeta_dr };

				// g' helper
				const auto dg_from = [&D, &denomG, &beta](const cplx& dbeta, const cplx& dD) noexcept -> cplx
					{
						return 2.0L * (D * dbeta - beta * dD) / denomG;
					};

				// Derivatives of g
				const cplx dg_dk{ dg_from(dbeta_dk, dD_dk) };
				const cplx dg_ds{ dg_from(dbeta_ds, dD_ds) };
				const cplx dg_dr{ dg_from(dbeta_dr, dD_dr) };

				// B' helper
				const auto dB_from = [&](const cplx& dbeta, const cplx& dD, long double dC) noexcept -> cplx
					{
						// Prefactor derivative wrt (beta, D, C=sigma)
						const cplx d_one_over_C2{ (-2.0L * dC) / sigma3 };
						const cplx dpref{ (dbeta - dD) * invSigma2 + betaMinusD * d_one_over_C2 };

						// Fraction derivative using cached inverses
						const cplx dN1{ (+deDT_dD) * (-dD) };
						const cplx dQ{ -(dg_from(dbeta, dD) * eDT) - g * (deDT_dD * dD) };
						const cplx dfrac{ (dN1 * Q - (1.0L - eDT) * dQ) * invQ2 };

						return dpref * fracB + betaMinusDinvSigma2 * dfrac;
					};

				// dC for each parameter (C := sigma)
				constexpr long double dC_dk{ 0.0L };
				constexpr long double dC_ds{ 1.0L };
				constexpr long double dC_dr{ 0.0L };

				// B partials
				const cplx dB_dk{ dB_from(dbeta_dk, dD_dk, dC_dk) };
				const cplx dB_ds{ dB_from(dbeta_ds, dD_ds, dC_ds) };
				const cplx dB_dr{ dB_from(dbeta_dr, dD_dr, dC_dr) };

				// A' pieces
				const cplx dK_dk{ theta * invSigma2 };
				const cplx dK_ds{ (-2.0L * kappaTheta) / sigma3 };
				constexpr cplx dK_dr{ 0.0L, 0.0L };

				// S' helper
				const auto dS_from = [&](const cplx& dbeta, const cplx& dD, const cplx& dg) noexcept -> cplx
					{
						const cplx dQ{ -(dg * eDT + g * (deDT_dD * dD)) };
						return (dbeta - dD) * T - 2.0L * (dQ * invQ + dg / R);
					};

				// S partials
				const cplx dS_dk{ dS_from(dbeta_dk, dD_dk, dg_dk) };
				const cplx dS_ds{ dS_from(dbeta_ds, dD_ds, dg_ds) };
				const cplx dS_dr{ dS_from(dbeta_dr, dD_dr, dg_dr) };

				// A partials
				const cplx dA_dk{ dK_dk * S + kFac * dS_dk };
				const cplx dA_ds{ dK_ds * S + kFac * dS_ds };
				const cplx dA_dr{ dK_dr * S + kFac * dS_dr };

				// Outputs (real parts)
				const cplx kernelPsi{ kernel * psi };
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
		const long double disc{ std::exp(-r * T) };
		const long double pref{ (F / std::numbers::pi_v<long double>) * std::exp(alpha * w) };
		const long double scale{ disc * pref };

		// Assemble price and gradients (unrolled)
		std::array<double, 6> out{};
		out[0] = static_cast<double>(disc * (R - pref * integrals[0]));
		out[1] = static_cast<double>(-scale * integrals[1]);
		out[2] = static_cast<double>(-scale * integrals[2]);
		out[3] = static_cast<double>(-scale * integrals[3]);
		out[4] = static_cast<double>(-scale * integrals[4]);
		out[5] = static_cast<double>(-scale * integrals[5]);
		return out;
	}

	template <std::size_t N>
	template <typename T>
	void Pricer<N>::setParams(const std::array<T, 5>& params) noexcept
	{
		params_ = Params
		{
			static_cast<long double>(params[0]), // kappa
			static_cast<long double>(params[1]), // theta
			static_cast<long double>(params[2]), // sigma
			static_cast<long double>(params[3]), // rho
			static_cast<long double>(params[4])  // v0 
		};
	}

	template <std::size_t N>
	long double Pricer<N>::getResidues(long double alpha,
		const long double F,
		const long double K) noexcept
	{
		if (alpha < -1.0) return F - K;
		else return 0.0L;
	}

	template <std::size_t N>
	long double Pricer<N>::getAlpha(long double w) const noexcept
	{
		if (w >= 0.0) return config_.alphaItm;
		else return config_.alphaOtm;
	}

	template <std::size_t N>
	long double Pricer<N>::getPhi(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double w) noexcept
	{
		if ((w * (rho - sigma * w / (v0 + kappa * theta * T))) >= 0.0L) return 0.0L;
		else return std::copysign(std::numbers::pi_v<long double> / 12.0L, w);
	}

	template <std::size_t N>
	cplx Pricer<N>::charFunction(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		const cplx& u) noexcept
	{
		// Define i 
		constexpr cplx i{ 0.0L, 1.0L };

		// u * ( u+ 1)
		cplx uu{ u * (u + i) };

		// sigma^2
		const long double sigma2{ sigma * sigma };

		// beta := kappa - i * sigma * rho * u
		const cplx beta{ kappa - i * sigma * rho * u };

		// D : = sqrt(beta^2 + sigma^2 * u * (u + i))
		const cplx D{ std::sqrt(beta * beta + sigma2 * uu) };

		// beta + D
		const cplx betaPlusD{ beta + D };

		// beta - D
		const cplx r
		{
			(std::real(beta * std::conj(D)) > 0.0L)
			? -sigma2 * uu / betaPlusD
			: beta - D
		};

		// y := expm1(−D * T) / (2 * D) with D≈0 fallback
		const cplx DT{ D * T };
		cplx y
		{
			(std::norm(D) > std::numeric_limits<long double>::epsilon() * (1.0L + std::abs(DT)))
			? math::expm1Complex(-DT) / (2.0L * D)
			: cplx(-T / 2.0L)
		};

		// r * y 
		const cplx ry{ r * y };

		// A := (κ * theta / sigma^2) * (r * T − 2 * log1p(−r * y))
		const cplx A{(kappa * theta / sigma2) * (r * T - 2.0L * math::log1pComplex<long double>(-ry))	};

		// B := u * (u + i) * y / (1 − r * y)
		const cplx B{ uu * y / (1.0L - ry) };

		// psi(u) := e^( A + v0 * B),
		return std::exp(A + v0 * B);
	}

	template <std::size_t N>
	CharFunData Pricer<N>::charFunctionCal(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		const cplx& u) noexcept
	{
		// i
		constexpr cplx i{ 0.0L, 1.0L };

		// u * (u + i)
		const cplx uu{ u * (u + i) };

		// sigma^2
		const long double sigma2{ sigma * sigma };

		// 1 / sigma^2
		const long double invSigma2{ 1.0L / sigma2 };

		// u * i
		const cplx ui{ u * i };

		// beta := kappa - sigma * rho * i*u
		const cplx beta{ kappa - sigma * rho * ui };

		// D := sqrt(beta^2 + sigma^2 * u(u+i))
		const cplx D{ std::sqrt(beta * beta + sigma2 * uu) };

		// beta + D
		const cplx betaPlusD{ beta + D };

		// beta - D  (stable branch)
		const cplx betaMinusD
		{
			(std::real(beta * std::conj(D)) > 0.0L)
			? -sigma2 * uu / betaPlusD
			: beta - D
		};

		// DT := D * T
		const cplx DT{ D * T };

		// y := expm1(-DT) / (2D)   with D≈0 fallback
		const cplx y
		{
			(std::norm(D) > std::numeric_limits<long double>::epsilon() *
							  (1.0L + std::abs(DT)))
			? math::expm1Complex(-DT) / (2.0L * D)
			: cplx(-T / 2.0L)
		};

		// r * y
		const cplx ry{ betaMinusD * y };

		// kappa * theta
		const long double kappaTheta{ kappa * theta };

		// (kappa * theta) / sigma^2
		const cplx kFac{ kappaTheta * invSigma2 };

		// A := (kappa*theta/sigma^2) * (r*T − 2 log(1 - r*y))
		const cplx A{ kFac * (betaMinusD * T - 2.0L * math::log1pComplex<long double>(-ry)) };

		// B := u(u+i)*y / (1 − r*y)
		const cplx B{ uu * y / (1.0L - ry) };

		// Rescued intermediates for gradient path

		// exp(-DT)
		const cplx eDT{ std::exp(-DT) };

		// g := (betaMinusD)/(betaPlusD)
		const cplx g{ betaMinusD / betaPlusD };

		// Q := 1 - g * eDT
		const cplx Q{ 1.0L - g * eDT };

		// 1/Q and 1/Q^2
		const cplx invQ{ 1.0L / Q };

		// R := 1 - g
		const cplx R{ 1.0L - g };

		// Return full CharFunData
		return
		{
			std::exp(A + v0 * B),                     // psi(u) := exp( A + v0 * B )
			A,                                          // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
			B,                                          // B := u(u+i)*y / (1 − r*y)
			beta,                                       // beta := kappa − sigma*rho*(i*u)
			D,                                          // D := sqrt( beta^2 + sigma^2*u(u+i) )
			DT,                                         // DT := D*T
			betaPlusD,                                  // beta + D
			betaMinusD,                                 // beta − D   (stable)
			ui,                                         // ui := u*i
			kFac,                                       // kFac := (kappa*theta)/sigma^2
			invSigma2,                                  // 1 / sigma^2
			kappaTheta,                                 // kappa * theta
			sigma2,                                     // sigma^2
			uu,                                         // uu := u(u+i)
			eDT,                                        // eDT := exp( −DT )
			betaMinusD / betaPlusD,                     // g := (beta−D)/(beta+D)
			Q,                                          // Q := 1 − g*eDT
			invQ,                                       // 1 / Q
			invQ * invQ,                                // 1 / Q^2
			R,                                          // R := 1 − g
			betaMinusD * T - 2.0L * std::log(Q / R),  // S := (beta−D)*T − 2*log(Q/R)
			(1.0L - eDT) * invQ,                        // fracB := (1 − eDT) / Q
			betaPlusD * betaPlusD,                      // denomG := (beta + D)^2
			betaMinusD * invSigma2                      // betaMinusD / sigma^2
		};
	}
}