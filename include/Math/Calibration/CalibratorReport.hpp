/**
* CalibratorReport.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/SVI/SVISlice.hpp"
#include "Core/VolSurface.hpp"
#include <vector>

namespace uv
{
	template<typename T>
	struct CalibratorReport
	{
		T params;               // Calibrated Report
		VolSurface volSurf;     // Calibrated volatility surface
	};
}