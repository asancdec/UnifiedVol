// SPDX-License-Identifier: Apache-2.0
/*
 * File:        HestonCalibrator.cpp
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


#include "Models/Heston/Calibrator.hpp"
#include "Utils/Aux/Errors.hpp"

#include <array>
#include <cmath>
#include <cstddef>

namespace uv::models::heston::calibrator::detail
{	
	std::array<double, 5> initGuess() noexcept
	{
		return
		{
			2.5,     // kappa
			0.09,    // theta
			0.60,    // sigma
			-0.75,   // rho
			0.09     // vo
		};
	}

	std::array<double, 5> lowerBounds() noexcept
	{
		return
		{
			0.001,   // kappa
			0.001,   // theta
			0.001,   // sigma
			-0.999,  // rho
			0.001    // vo
		};

	}

	std::array<double, 5> upperBounds() noexcept
	{
		return
		{
			10.0,   // kappa
			0.5,    // theta
			10.0,   // sigma
			0.999,  // rho
			0.5     // vo
		};
	}
} // namespace uv::models::heston_calibrator::detail

