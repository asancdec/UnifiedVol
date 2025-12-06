/**
* LocalVol.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"
#include "Utils/Types/Types.hpp"
#include "Models/SVI/SVISlice.hpp"

namespace uv::localvol
{	
	// Build local volatility surface from SVI parameters
	VolSurface buildSurface(const VolSurface& volSurface, 
		const Vector<models::svi::SVISlice>& sviSlices);

	namespace detail
	{	
		// Calculate total variance first derivative w.r.t. time using piecewise derivatives
		Matrix<double> dwdT(const Vector<double>& tenors,
			const Matrix<double>& totVarMatrix) noexcept;
	}

} // uv::localvol
