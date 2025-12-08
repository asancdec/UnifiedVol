// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Config.hpp
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

#include "Utils/Types.hpp"
#include <array>    
#include <string_view>
#include <cstddef>

namespace uv::math::opt::nlopt
{
	template <std::size_t N>
	struct Config
	{
		double eps;                                       // Small epsilon used for inequality tolerances
		double tol;                                       // Constraint tolerance passed to NLopt (|c(x)| ≤ tol considered satisfied)
		double ftolRel;                                   // Relative stopping criterion for objective improvement (NLopt ftol_rel)
		unsigned int maxEval;                             // Maximum number of function evaluations allowed during optimization
		std::array<std::string_view, N> paramNames;       // Parameter names
	};
}
