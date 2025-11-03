/**
* HestonPricer.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Errors/Errors.hpp"
#include "Math/MathFunctions/MathFunctions.hpp"

#include <utility>   
#include <complex>  
#include <cmath>     
#include <numbers>  
#include <limits>    

namespace uv
{
	template <::std::size_t N>
	HestonPricer<N>::HestonPricer(::std::shared_ptr<const TanHSinH<N>> quad,
		const HestonConfig& config) :
		quad_(::std::move(quad)),
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

	template <::std::size_t N>
	double HestonPricer<N>::callPrice(long double kappa,
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
		const long double w{ ::std::log(F / K) };

		// Determine alpha
		const long double alpha(getAlpha(w));

		// Calculate residues
		const long double R{ HestonPricer<N>::getResidues(alpha, F, K) };

		// Determine phi
		const long double phi{ HestonPricer<N>::getPhi(kappa, theta, sigma, rho, v0, T, w) };

		// Precomputations
		constexpr ::std::complex<long double> i(0.0L, 1.0L);
		const long double tanPhi{ ::std::tan(phi) };
		const ::std::complex<long double> onePlusITanPhi{ 1.0L + i * tanPhi };
		const ::std::complex<long double> c{ (i - tanPhi) * w };

		// Define the integrand
		auto integrand = [=](long double x) noexcept -> long double
			{
				// Calculate h(x)
				const ::std::complex<long double> h{ -i * alpha + x * onePlusITanPhi };

				// h - i
				const ::std::complex<long double> hMinusI{ h - i };

				// Evaluate characteristic function at h(x) - i
				const ::std::complex<long double> charFuncVal{ HestonPricer<N>::charFunction(kappa, theta, sigma, rho, v0, T, hMinusI) };

				// Calculate Q(h(x))
				const ::std::complex<long double> Q{ charFuncVal / (hMinusI * h) };

				// Calculate and return integrand
				return ::std::real
				(
					::std::exp(x * c) * Q * onePlusITanPhi
				);
			};

		// Evaluate integrand
		const long double integral{ quad_->integrateZeroToInf(integrand) };

		// Calculate and return Call price
		return static_cast<double>(::std::exp(-r * T) * (R - F / ::std::numbers::pi_v<long double> *::std::exp(alpha * w) * integral));
	}

	template <::std::size_t N>
	double HestonPricer<N>::callPrice(long double T,
		long double F,
		long double r,
		long double K) const
	{
		// Will return error if class parameters are not set
		UV_REQUIRE
		(
			(params_.has_value()),
			ErrorCode::InvalidArgument,
			"HestonPricer::callPrice: Heston parameters not set. Call setHestonParams(...) first."
		);

		// Dereference optional
		const HestonParams& params{ *params_ };

		// Pass class instance parameters into the generic pricing function
		return callPrice(params.kappa,
			params.theta,
			params.sigma,
			params.rho,
			params.v0,
			T, F, r, K);
	}


	template <::std::size_t N>
	long double HestonPricer<N>::getResidues(long double alpha,
		const long double F,
		const long double K) noexcept
	{
		if (alpha < -1.0) return F - K;
		else return 0.0L;
	}

	template <::std::size_t N>
	long double HestonPricer<N>::getAlpha(long double w) const noexcept
	{
		if (w >= 0.0) return config_.alphaItm;
		else return config_.alphaOtm;
	}

	template <::std::size_t N>
	long double HestonPricer<N>::getPhi(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double w) noexcept
	{
		if ((w * (rho - sigma * w / (v0 + kappa * theta * T))) >= 0.0L) return 0.0L;
		else return ::std::copysign(::std::numbers::pi_v<long double> / 12.0L, w);
	}

	template <::std::size_t N>
	::std::complex<long double> HestonPricer<N>::charFunction(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		::std::complex<long double> u) noexcept
	{
		// Define i 
		constexpr ::std::complex<long double> i(0.0L, 1.0L);

		// u * ( u+ 1)
		::std::complex<long double> uu{ u * (u + ::std::complex<long double>(0.0L, 1.0L)) };

		// sigma^2
		const long double sigma2{ sigma * sigma };

		// beta := kappa - i * sigma * rho * u
		const ::std::complex<long double> beta{ kappa - ::std::complex<long double>(0.0L, 1.0L) * sigma * rho * u };

		// D : = sqrt(beta^2 + sigma^2 * u * (u + i))
		const ::std::complex<long double> D{ ::std::sqrt(beta * beta + sigma2 * uu) };

		// beta + D
		const ::std::complex<long double> betaPlusD{ beta + D };

		// beta - D
		const ::std::complex<long double> r
		{
			(::std::real(beta * ::std::conj(D)) > 0.0L)
			? -sigma2 * uu / betaPlusD
			: beta - D
		};

		// y := expm1(−D * T) / (2 * D) with D≈0 fallback
		const ::std::complex<long double> DT{ D * T };
		::std::complex<long double> y
		{
			(::std::abs(D) > ::std::numeric_limits<long double>::epsilon() * (1.0L + ::std::abs(DT)))
			? expm1Complex(-DT) / (2.0L * D)
			: ::std::complex<long double>(-T / 2.0L)
		};

		// r * y 
		const ::std::complex<long double> ry{ r * y };

		// A := (κ * theta / sigma^2) * (r * T − 2 * log1p(−r * y))
		const ::std::complex<long double> A
		{
			(kappa * theta / sigma2) * (r * T - 2.0L * log1pComplex<long double>(-ry))
		};

		// B := u * (u + i) * y / (1 − r * y)
		const ::std::complex<long double> B{ uu * y / (1.0L - ry) };

		// psi(u) := e^( A + v0 * B)
		return ::std::exp(A + v0 * B);
	}

	template <::std::size_t N>
	template <typename T>
	void HestonPricer<N>::setHestonParams(const ::std::array<T, 5>& params) noexcept
	{	
		// Set pricer parameters
		params_ = HestonParams
		{
			static_cast<long double>(params[0]), // kappa
			static_cast<long double>(params[1]), // theta
			static_cast<long double>(params[2]), // sigma
			static_cast<long double>(params[3]), // rho
			static_cast<long double>(params[4])  // v0 
		};
	}

}