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
#include <memory>
#include <optional>    
#include <cstddef>  
#include <span>

namespace uv
{	
	template <::std::size_t N>
	class HestonPricer
	{
	private:

		//--------------------------------------------------------------------------
		// Private member variables
		//--------------------------------------------------------------------------
		::std::optional<HestonParams> params_;
		::std::shared_ptr<const TanHSinH<N>> quad_;
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
		explicit HestonPricer(::std::shared_ptr<const TanHSinH<N>> quad,
			const HestonConfig& config);

		//--------------------------------------------------------------------------
		// Pricing
		//--------------------------------------------------------------------------
		// Overload 1: using user-defined parameters (during the calibration for example) 
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

		// Overload 2: using class parameters (typically already calibrated) 
		// Andersen & Lake Implementation
		double callPrice(long double T,
			long double F,
			long double r,
			long double K) const;

		//--------------------------------------------------------------------------
		// Setters
		//--------------------------------------------------------------------------
		template <typename T>
		void setHestonParams(const ::std::array<T, 5>& params) noexcept;
	};
}

#include "HestonPricer.inl"