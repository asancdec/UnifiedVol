// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */


#pragma once

#include "Core/VolSurface.hpp"
#include "Utils/Types.hpp"
#include "Models/SVI/Params.hpp"

#include <span>

namespace uv::models::localvol
{	
	/**
	 * @brief Build a local volatility surface using SVI slices parameters.
	 *
	 * This function takes an existing total variance volatility surface and a set
	 * of SVI parameter slices (one per tenor), and builds the local total
	 * variance surface using Gatheral–style local volatility construction.
	 *
	 * @param volSurface
	 *     Input total variance surface
	 *
	 * @param paramsSVI
	 *     Vector of SVI parameter slices. Must have same tenors as volSurface
	 *
	 * @return core::VolSurface
	 *     A new surface with identical tenor/strike structure to the input
	 *     surface, but with each slice's total variance replaced by the local
	 *     total variance computed via SVI..
	 */
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
		const Vector<models::svi::Params>& paramsSVI);

	/**
	 * @brief Price a European call option surface using a local volatility model.
	 *
	 * This function constructs a local volatility pricer with the supplied market data
	 * and grid configuration, and computes call option prices for the full set of
	 * strikes and maturities defined in the input volatility surface.
	 *
	 * @param localVolSurface  Calibrated local volatility surface (tenors, strikes, vol matrix).
	 * @param marketData       Spot, rate, and dividend yield used for pricing.
	 * @param NT               Number of time grid points for the PDE solver.
	 * @param NS               Number of spot grid points for the PDE solver.
	 * @param X                Spot upper-bound multiplier (domain scaling parameter).
	 *
	 * @return Vector<Real>    Matrix-flattened call prices consistent with the given surface.
	 */
	Vector<Real> priceCall(const core::VolSurface& localVolSurface,
		const core::MarketData& marketData,
		const std::size_t NT,
		const std::size_t NS,
		const unsigned int X = 3);

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
		 * @param paramsSVI
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
            const Vector<svi::Params>& paramsSVI);
	}

} // uv::models::localvol
