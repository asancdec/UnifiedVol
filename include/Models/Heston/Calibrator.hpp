/**
* Calibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Models/Heston/Pricer.hpp"
#include "Models/Heston/Params.hpp"
#include "Calibration/Ceres/Calibrator.hpp"
#include "Core/VolSurface.hpp"

#include <array>
#include <cstddef>  

namespace uv::models::heston::calibrator
{
	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Main calibration function
	template <std::size_t N, typename Policy>
	static core::VolSurface calibrate(const core::VolSurface& mktVolSurf,
		Pricer<N>& pricer,
		cal::ceres::Calibrator<5, Policy>& calibrator);

	namespace detail
	{
		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------

		// Residue functor per call price
		template <std::size_t N>
		struct PriceResidualJac;

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------

		// Initial guess
		std::array<double, 5> initGuess() noexcept;

		// Lower bounds
		std::array<double, 5> lowerBounds() noexcept;

		// Upper bounds
		std::array<double, 5> upperBounds() noexcept;
	} //  namespace detail
}

#include "Calibrator.inl"