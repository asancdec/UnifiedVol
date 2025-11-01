/**
* HestonPricer.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/Heston/HestonParams.hpp"
#include "Models/Heston/HestonConfig.hpp"
#include "Math/Quadrature/TanHSinH.hpp"
#include "Core/VolSurface.hpp"
#include <complex>
#include <cstdint> 
#include <memory>
#include <optional>

namespace uv
{	
	class HestonPricer
	{
	private:

		//--------------------------------------------------------------------------
		// Private member variables
		//--------------------------------------------------------------------------
		::std::optional<HestonParams> params_;
		::std::shared_ptr<const TanHSinH> quad_;
		const HestonConfig config_;

		//--------------------------------------------------------------------------
		// Math
		//--------------------------------------------------------------------------
		// Calculate residues arising from the contour shift
		static long double getResidues(long double alpha,
			long double F,
			long double K) noexcept;

		// Get ITM or OTM damping parameter
		long double getAlpha(long double w) const noexcept;

		// Determine contour shift angle
		static long double getPhi(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			long double w) noexcept;

		// No branch cut characteristic function 
		// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
		static ::std::complex<long double> charFunction(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			::std::complex<long double> u) noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		HestonPricer() = delete;
		explicit HestonPricer(::std::shared_ptr<const TanHSinH> quad,
			const HestonConfig& config);

		//--------------------------------------------------------------------------
		// Pricing
		//--------------------------------------------------------------------------
		// European Call Price
		double callPrice(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			long double F,
			long double r,
			long double K) const noexcept;
	};
}