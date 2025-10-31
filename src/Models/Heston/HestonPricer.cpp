/**
* Heston.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/Heston/HestonPricer.hpp"
#include "Errors/Errors.hpp"
#include "Math/MathFunctions/MathFunctions.hpp"
#include <cmath>
#include <numbers>
#include <iostream>


namespace uv
{
	HestonPricer::HestonPricer(std::shared_ptr<const TanHSinH> quad,
		const HestonConfig& config) :
		quad_(std::move(quad)),
		config_(config)
	{
		// Alpha checks
		UV_REQUIRE(
			(config_.alphaItm <= -1.0 - config_.eps) && (config_.alphaOtm >= config_.eps),
			ErrorCode::InvalidArgument,
			"Alpha must be outside [-1, 0] range");

		std::cout << "Heston started\n";
	}

	double HestonPricer::callPrice(long double kappa,
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
		const long double R{ getResidues(alpha, F, K) };

		// Determine phi
		const long double phi{ getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr std::complex<long double> i(0.0L, 1.0L);
		const long double tanPhi{ std::tan(phi) };
		const std::complex<long double> onePlusITanPhi{ 1.0L + i * tanPhi };
		const std::complex<long double> c{ (i - tanPhi) * w };

		// Define the integrand
		auto integrand = [=](long double x) noexcept -> long double
			{
				// Calculate h(x)
				const std::complex<long double> h{ -i * alpha + x * onePlusITanPhi };

				// h - i
				const std::complex<long double> hMinusI{ h - i };

				// Evaluate characteristic function at h(x) - i
				const std::complex<long double> charFuncVal{ charFunction(kappa, theta, sigma, rho, v0, T, hMinusI) };

				// Calculate Q(h(x))
				const std::complex<long double> Q{ charFuncVal / (hMinusI * h) };

				// Calculate and return integrand
				return std::real
				(
					std::exp(x * c) * Q * onePlusITanPhi
				);
			};

		// Evaluate integrand
		const long double integral{ quad_->integrateZeroToInf(integrand) };

		// Calculate and return Call price
		return static_cast<double>(std::exp(-r * T) * (R - F / std::numbers::pi_v<long double> *std::exp(alpha * w) * integral));
	}

	long double HestonPricer::getResidues(long double alpha,
		const long double F,
		const long double K) noexcept
	{
		if (alpha < -1.0) return F - K;
		else return 0.0L;
	}

	long double HestonPricer::getAlpha(long double w) const noexcept
	{
		if (w >= 0.0) return config_.alphaItm;
		else return config_.alphaOtm;
	}

	long double HestonPricer::getPhi(long double kappa,
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

	std::complex<long double> HestonPricer::charFunction(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		std::complex<long double> u) noexcept
	{
		// Define i 
		constexpr std::complex<long double> i(0.0L, 1.0L);

		// u * ( u+ 1)
		std::complex<long double> uu{ u * (u + i) };

		// sigma^2
		const long double sigmaSquared{ sigma * sigma };

		// beta := kappa - i * sigma * rho * u
		const std::complex<long double> beta{ kappa - i * sigma * rho * u };

		// D : = sqrt(beta^2 + sigma^2 * u * (u + i))
		const std::complex<long double> D{ std::sqrt(beta * beta + sigmaSquared * uu) };

		// B + D
		const std::complex<long double> betaPlusD{ beta + D };

		// B - D
		const std::complex<long double> betaMinusD{ -sigmaSquared * uu / betaPlusD };

		// G := (beta - D) / (beta + D)
		const std::complex<long double> G{ betaMinusD / betaPlusD };

		// e^(-D * T)
		const std::complex<long double> expMinusDT{ std::exp(-D * T) };

		// 1 - G * e^(-D*T)
		const std::complex<long double> oneMinusGExpMinusDt{ 1.0L - G * expMinusDT };

		// A : = k * theta / sigma^2 * ((beta - D) * T - 2 * ln((1 - G * e^(-D *  T)) / (1 - G)))
		const std::complex<long double> A{ kappa * theta / sigmaSquared * (betaMinusD * T - 2.0L
			* (log1pComplex<long double>(-G * expMinusDT - log1pComplex<long double>(-G)))) };

		// B: = (beta - D) / sigma^2 * ((1 - e^(-D * T)) / (1 - G * e^(-D * T)))
		const std::complex<long double>B{ betaMinusD / sigmaSquared * (1.0L - expMinusDT) / oneMinusGExpMinusDt };

		// psi(u) := e^( A + v0 * B)
		return std::exp(A + v0 * B);
	}
}
