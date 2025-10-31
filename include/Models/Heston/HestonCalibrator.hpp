/**
* HestonCalibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/Heston/HestonPricer.hpp"
#include "Models/Heston/HestonParams.hpp"
#include "Math/Calibration/CalibratorReport.hpp"
#include "Core/VolSurface.hpp"
#include <array>

namespace uv
{
	using HestonReport = CalibratorReport<HestonParams>;

	class HestonCalibrator
	{
	private:

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------

		// Initial guess
		static std::array<double, 5> initGuess() noexcept;

		// Lower bounds
		static std::array<double, 5> lowerBounds() noexcept;

		// Upper bounds
		static std::array<double, 5> upperBounds() noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		HestonCalibrator() = default;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		static HestonReport calibrate(const VolSurface mktVolSurf,
			const HestonPricer& pricer);
	};
}
