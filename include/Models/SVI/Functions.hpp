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
#include "Math/Optimization/NLopt/Optimizer.hpp"
#include "Models/SVI/Params.hpp"
#include "Utils/Types.hpp"

#include <nlopt.hpp> 
#include <concepts>
#include <array>
#include <vector>
#include <tuple>

namespace opt = uv::math::opt;

namespace uv::models::svi
{
	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Main calibration function
	template <::nlopt::algorithm Algo>
	std::tuple<std::vector<Params>, core::VolSurface> calibrate(const core::VolSurface& mktVolSurf,
		const opt::nlopt::Optimizer<4, Algo>& prototype);

	//// g(k) standard function
	//double gk(Real a, Real b, Real rho, Real m, Real sigma, Real k) noexcept;

	namespace detail
	{
		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------

		struct CalendarContexts;
		struct ObjectiveContexts;
		struct GkCache;
		struct ConvexityContexts;

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------

		// Initial guess
		std::array<double, 4> initGuess() noexcept;

		// Lower bounds
		std::array<double, 4> lowerBounds(const double logKFMin) noexcept;

		// Upper bounds
		std::array<double, 4> upperBounds(const double logKFMax) noexcept;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------

		// Define the minimum total variance constraint: wMin ≥ 0.0
		template <::nlopt::algorithm Algo>
		void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

		// Add Roger Lee left wing and right wing min slope constraints
		template <::nlopt::algorithm Algo>
		void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

		// Add Roger Lee left wing and right wing max slope constraints
		template <::nlopt::algorithm Algo>
		void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

		// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
		template <::nlopt::algorithm Algo>
		void addCalendarConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer, Vector<CalendarContexts>& ctx);

		// Add no butterfly arbitrage constraints g(k) ≥ 0.0
		template <::nlopt::algorithm Algo>
		void addConvexityConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer, ConvexityContexts& ctx);

		// Add objective function with analytical gradient
		template <::nlopt::algorithm Algo>
		void setMinObjective(opt::nlopt::Optimizer<4, Algo>& optimizer, const ObjectiveContexts& ctx);

		//--------------------------------------------------------------------------
		// Math functions
		//--------------------------------------------------------------------------

		// SVI total variance for a given log-forward moneyness
		double calculateWk(const double a,
			const double b,
			const double rho,
			const double m,
			const double sigma,
			const double k) noexcept;

		// a param from the atmWT
		double aParam(const double atmWK,
			const double b,
			const double rho,
			const double m,
			const double sigma) noexcept;

		// g(k) optimized function for calibration
		double gk(const GkCache& p) noexcept;

		// ∇g(k)
		std::array<double, 4> gkGrad(const double b,
			const double rho,
			const double m,
			const double sigma,
			const double k,
			const GkCache& p) noexcept;

		// Make a wK slice
		template <std::floating_point T>
		Vector<T> makewKSlice(const Vector<T>& kSlice,
			T a, T b, T rho, T m, T sigma) noexcept;

	} // namespace detail
} 

#include "Functions.inl"