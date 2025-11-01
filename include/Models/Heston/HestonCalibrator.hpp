/**
* HestonCalibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/Heston/HestonPricer.hpp"
#include "Models/Heston/HestonParams.hpp"
#include "Math/Calibration/CalibratorReport.hpp"
#include "Math/Calibration/Ceres/CalibratorCeres.hpp"
#include "Core/VolSurface.hpp"
#include <array>
#include <memory>

namespace uv
{
	using HestonReport = CalibratorReport<HestonParams>;

	class HestonCalibrator
	{
	private:

		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------
		// Residue functor per call price
		// Uses numeric differentiation
		struct PriceResidualND;

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
		static void calibrate(const VolSurface& mktVolSurf,
			::std::shared_ptr<const HestonPricer> pricer,
			CalibratorCeres<5>& calibrator);
	};
}
