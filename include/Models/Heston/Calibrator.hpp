// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Calibrator.hpp
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

#include "Models/Heston/Pricer.hpp"
#include "Models/Heston/Params.hpp"
#include "Math/Optimization/Ceres/Optimizer.hpp"
#include "Core/VolSurface.hpp"
#include "Core/Matrix.hpp"

#include <array>
#include <cstddef>

namespace uv::models::heston::calibrator
{
	namespace opt = uv::math::opt;

	/**
	 * @brief Calibrate Heston parameters to a matrix of market call prices using Ceres.
	 *
	 * The calibration minimizes per-instrument pricing residuals:
	 *   residual(T_i, K_j) = modelCall(T_i, F_i, r_i, K_j) - marketCall(T_i, K_j)
	 *
	 * The optimizer is configured with an initial guess and bounds, then one analytic
	 * residual block is added for each (tenor, strike) point.
	 *
	 * @tparam N Quadrature setting used by the pricer backend.
	 * @tparam Policy Optimizer policy used by the Ceres wrapper (loss/solver settings).
	 *
	 * @param tenors   Tenor grid (size = numTenors), strictly increasing, > 0.
	 * @param strikes  Strike grid (size = numStrikes), strictly increasing, > 0.
	 * @param forwards Forward values per tenor (size = numTenors), > 0.
	 * @param rates    Discounting rates per tenor (size = numTenors).
	 * @param callM    Market call price matrix [numTenors x numStrikes].
	 * @param pricer   Heston pricer used to evaluate prices and analytic gradients.
	 * @param optimizer Ceres optimizer wrapper for 5 Heston parameters.
	 *
	 * @return Calibrated Heston parameters (kappa, theta, sigma, rho, v0).
	 *
	 * @throws (via UV_REQUIRE) if input dimensions are inconsistent or invalid.
	 */
	template <std::size_t N, typename Policy>
	Params calibrate(const Vector<Real>& tenors,
		const Vector<Real>& strikes,
		const Vector<Real>& forwards,
		const Vector<Real>& rates,
		const core::Matrix<Real>& callM,
		Pricer<N>& pricer,
		opt::ceres::Optimizer<5, Policy>& optimizer);

	/**
	 * @brief Build a new volatility surface by repricing all options with the calibrated Heston pricer.
	 *
	 * This function:
	 *  1) reads (tenors, strikes, forwards, rates) from the input surface,
	 *  2) computes a model call price matrix using pricer.callPrice(T, F, r, K),
	 *  3) returns a copy of the original surface with call prices replaced.
	 *
	 * @tparam N Quadrature setting used by the pricer backend.
	 *
	 * @param volSurface Input surface providing grids and market metadata.
	 * @param pricer     Heston pricer (must have parameters set, or be able to price).
	 *
	 * @return Copy of volSurface with call prices set to the Heston model prices.
	*/
	template <std::size_t N>
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
		const Pricer<N>& pricer);

	namespace detail
	{
		/**
		 * @brief Validate calibration inputs for size, monotonicity, and numeric sanity.
		 *
		 * Checks:
		 *  - non-empty tenors/strikes/forwards/rates/callM,
		 *  - size consistency: forwards/rates match tenors, callM rows match tenors,
		 *    callM columns match strikes,
		 *  - tenors strictly increasing and > 0,
		 *  - strikes strictly increasing and > 0,
		 *  - forwards > 0,
		 *  - callM rectangular and finite.
		 *
		 * @param tenors   Tenor grid.
		 * @param strikes  Strike grid.
		 * @param forwards Forwards per tenor.
		 * @param rates    Rates per tenor.
		 * @param callM    Market call matrix.
		 *
		 * @throws (via UV_REQUIRE) if any condition fails.
		 */
		void validateInputs(const Vector<Real>& tenors,
			const Vector<Real>& strikes,
			const Vector<Real>& forwards,
			const Vector<Real>& rates,
			const core::Matrix<Real>& callM);

		/**
		 * @brief Ceres analytic residual for a single call price point.
		 *
		 * Computes:
		 *  - residual = modelPrice - marketPrice
		 * and (optionally) the Jacobian w.r.t. the 5 Heston parameters:
		 *  [kappa, theta, sigma, rho, v0].
		 *
		 * @tparam N Quadrature setting used by the pricer backend.
		 */
		template <std::size_t N>
		struct PriceResidualJac;

		/**
		 * @brief Default initial guess for Heston parameters used by the calibrator.
		 * @return Array in order: [kappa, theta, sigma, rho, v0].
		 */
		std::array<double, 5> initGuess() noexcept;

		/**
		 * @brief Default lower bounds for Heston parameters used by the calibrator.
		 * @return Array in order: [kappa, theta, sigma, rho, v0].
		 */
		std::array<double, 5> lowerBounds() noexcept;

		/**
		 * @brief Default upper bounds for Heston parameters used by the calibrator.
		 * @return Array in order: [kappa, theta, sigma, rho, v0].
		 */
		std::array<double, 5> upperBounds() noexcept;

	} //  namespace detail
} // namespace  uv::models::heston::calibrator

#include "Calibrator.inl"