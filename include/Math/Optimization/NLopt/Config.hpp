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

#include <array>
#include <cstddef>
#include <string_view>

namespace uv::math::opt::nlopt
{
/**
 * @brief Configuration parameters for an NLopt optimizer.
 *
 * Holds numerical tolerances and metadata used to configure NLopt-based
 * optimization routines.
 *
 * This structure does not own any optimizer state; it simply groups
 * commonly used stopping criteria and constraint tolerances in one place.
 *
 * @tparam N Number of optimization parameters.
 */
template <std::size_t N> struct Config
{
    double eps;           ///< Small epsilon used in inequality constraints
    double tol;           ///< Constraint satisfaction tolerance
    double ftolRel;       ///< Relative tolerance for objective improvement
    unsigned int maxEval; ///< Maximum number of objective evaluations
    std::array<std::string_view, N>
        paramNames; ///< Human-readable parameter names (for logging/debugging)
};
} // namespace uv::math::opt::nlopt
