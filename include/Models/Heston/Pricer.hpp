/**
* Pricer.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Models/Heston/Params.hpp"
#include "Models/Heston/Config.hpp"
#include "Models/Heston/CharFunData.hpp"
#include "Math/Quadratures/TanHSinH.hpp"
#include "Core/VolSurface.hpp"

#include <complex>     
#include <memory>
#include <array>
#include <optional>    
#include <cstddef>  
#include <tuple>


namespace uv::models::heston
{	

	template <std::size_t N>
	class Pricer
	{
	private:

		//--------------------------------------------------------------------------
		// Private member variables
		//--------------------------------------------------------------------------
		std::optional<Params> params_;
		std::shared_ptr<const math::TanHSinH<N>> quad_;
		const Config config_;

		//--------------------------------------------------------------------------
		// Pricing
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

		// Returns only price
		// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
		static cplx charFunction(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			const cplx& u) noexcept;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		// Returns a struct of precalculated variables for efficient gradient computation
		// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
		static CharFunData charFunctionCal(long double  kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			const cplx& u) noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		Pricer() = delete;
		explicit Pricer(std::shared_ptr<const math::TanHSinH<N>> quad,
			const Config& config);

		//--------------------------------------------------------------------------
		// Pricing
		//--------------------------------------------------------------------------
		// Overload 1: using user-defined parameters
		// Andersen & Lake Implementation
		double callPrice(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			long double F,
			long double r,
			long double K) const noexcept;

		// Overload 2: using class instance parameters
		// Andersen & Lake Implementation
		double callPrice(long double T,
			long double F,
			long double r,
			long double K) const;

		// Calculate price and parameter gradient for the calibration
		std::array<double, 6> callPriceWithGradient(long double kappa,
			long double theta,
			long double sigma,
			long double rho,
			long double v0,
			long double T,
			long double F,
			long double r,
			long double K) const noexcept;

		//--------------------------------------------------------------------------
		// Setters
		//--------------------------------------------------------------------------
		template <typename T>
		void setParams(const std::array<T, 5>& params) noexcept;
	};
}

#include "Pricer.inl"