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

#include <array>
#include <cstddef>  

namespace uv::models::heston::calibrator
{

	namespace opt = uv::math::opt;

	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Main calibration function
	template <std::size_t N, typename Policy>
	core::VolSurface calibrate(const core::VolSurface& mktVolSurf,
		Pricer<N>& pricer,
		opt::ceres::Optimizer<5, Policy>& optimizer);

	namespace detail
	{
		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------

		// Residue functor per call price
		template <std::size_t N>
		struct PriceResidualJac;

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------

		// Initial guess
		std::array<double, 5> initGuess() noexcept;

		// Lower bounds
		std::array<double, 5> lowerBounds() noexcept;

		// Upper bounds
		std::array<double, 5> upperBounds() noexcept;
	} //  namespace detail
} // namespace  uv::models::heston::calibrator

#include "Calibrator.inl"