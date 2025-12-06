/**
* LocalVol.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"
#include "Utils/Types.hpp"
#include "Models/SVI/SVISlice.hpp"

#include <span>

namespace uv::models::localvol
{	
	/**
	 * @brief Build a local volatility surface using SVI slices.
	 *
	 * This function takes an existing total variance volatility surface and a set
	 * of SVI parameter slices (one per tenor), and builds the local total
	 * variance surface using Gatheral–style local volatility construction.
	 *
	 * @param volSurface
	 *     Input total variance surface
	 *
	 * @param sviSlices
	 *     Vector of SVI parameter slices. Must have same tenors as volSurface
	 *
	 * @return core::VolSurface
	 *     A new surface with identical tenor/strike structure to the input
	 *     surface, but with each slice's total variance replaced by the local
	 *     total variance computed via SVI..
	 */
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
		const Vector<models::svi::SVISlice>& sviSlices);

	namespace detail
	{	
		/**
		 * @brief Compute local total variance from SVI parameters and total variance.
		 *
		 * Implements the Gatheral formula for local volatility in terms of local total
		 * variance:
		 *
		 *        w_local(T, k) = (∂w/∂T)(T, k) * T / g(k),
		 *
		 * where:
		 *  - w(T, k) is total implied variance,
		 *  - g(k) is the SVI convexity denominator defined by `svi::gk(...)`,
		 *  - ∂w/∂T(T, k) is computed using a monotone cubic Hermite spline (PCHIP)
		 *    interpolation along the time dimension.
		 *
		 * @param tenors
		 *     Vector of tenor values of size N.
		 *
		 * @param logFM
		 *     N×M matrix of log-moneyness values log(F/K). Dimensions must match
		 *     `totVar`.
		 *
		 * @param totVar
		 *     N×M matrix of total implied variance w(T, k).
		 *
		 * @param sviSlices
		 *     Vector of SVI parameter sets, one per tenor. Must have size N, and
		 *     each slice provides (a, b, rho, m, sigma, T).
		 *
		 * @return Matrix<Real>
		 *     N×M matrix of local total variance values w_local(T, k), same structure
		 *     as the input matrices.
		 *
		 * @details
		 *     No extrapolation is performed beyond the provided tenor grid.
		 */
        Matrix<Real> localTotVar(const Vector<Real>& tenors,
            const Matrix<Real>& logFM,
            const Matrix<Real>& totVarMatrix,
            const Vector<svi::SVISlice>& sviSlices);
	}

} // uv::models::localvol
