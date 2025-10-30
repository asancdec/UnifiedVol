/**
* SVIReport.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/SVI/SVISlice.hpp"
#include "Core/VolSurface.hpp"
#include <vector>

struct SVIReport
{
	std::vector<SVISlice> sviSlices;   // Vector with calibrated SVI params of every slice
	VolSurface volSurf;				   // Calibrated volatility surface
};
