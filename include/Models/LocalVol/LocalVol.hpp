/**
* LocalVol.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVISlice.hpp"

namespace uv::local_vol
{	
	// Build local volatility surface from SVI parameters
	VolSurface build(const VolSurface& volSurface, const std::vector<SVISlice>& sviSlices);

	namespace detail
	{	
		// Calculate total variance first derivative w.r.t. time using piecewise derivatives
		std::vector<std::vector<double>> dwdT(const std::vector<double>& tenors,
			const std::vector<std::vector<double>>& totVarMatrix) noexcept;
	}

} // uv::local_vol
