/**
* HestonCalibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Models/Heston/HestonPricer/HestonPricer.hpp"
#include "Models/Heston/HestonParams.hpp"
#include "Math/Calibration/Ceres/CalibratorCeres.hpp"
#include "Core/VolSurface.hpp"

#include <array>
#include <cstddef>  

namespace uv
{
	class HestonCalibrator
	{
	private:

		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------
		// Residue functor per call price
		template <::std::size_t N>
		struct PriceResidualJac;

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------
		// Initial guess
		static ::std::array<double, 5> initGuess() noexcept;

		// Lower bounds
		static ::std::array<double, 5> lowerBounds() noexcept;

		// Upper bounds
		static ::std::array<double, 5> upperBounds() noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		HestonCalibrator() = default;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		template <::std::size_t N, typename Policy>
		static VolSurface calibrate(const VolSurface& mktVolSurf,
			HestonPricer<N>& pricer,
			CalibratorCeres<5, Policy>& calibrator);
	};
}

#include "HestonCalibrator.inl"