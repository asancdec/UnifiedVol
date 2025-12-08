// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Optimizer.hpp
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

#include "Math/Optimization/Ceres/Config.hpp"
#include "Policy.hpp"

#include <ceres/ceres.h>
#include <array>

namespace uv::math::opt::ceres
{
	template <std::size_t N, typename Policy = Policy<>>
	class Optimizer
	{
	private:

		//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------
		Config<N> config_;                   // Optimizer configuration
		std::array<double, N> lowerBounds_;  // Lower parameter bounds
		std::array<double, N> upperBounds_;  // Upper parameter bounds
		std::array<double, N> x_;            // Parameter block
		::ceres::Problem problem_;           // Ceres problem instance

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		Optimizer() = delete;
		explicit Optimizer(const Config<N>& config);

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------	
		// Set initial guess and bounds
		void setGuessBounds(const std::array<double, N>& initGuess,
			const std::array<double, N>& lowerBounds,
			const std::array<double, N>& upperBounds) noexcept;
		
		// Set differentiation cost function
		void addAnalyticResidual(std::unique_ptr<::ceres::CostFunction> cf) noexcept;

		// Solve the problem
		std::array<double, N> optimize();
	};
}

#include "Optimizer.inl"