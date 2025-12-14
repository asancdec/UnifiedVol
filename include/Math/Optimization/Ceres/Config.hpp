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
#include <string_view>
#include <cstddef>

namespace uv::math::opt::ceres
{
    /**
     * @brief Configuration for Ceres-based optimisation.
     *
     * Groups solver tolerances, iteration limits, and logging options used by
     * the Ceres optimiser wrapper.
     *
     * @tparam N Number of optimisation parameters.
     */
    template <std::size_t N>
    struct Config
    {
        unsigned maxEval;                               // Maximum number of iterations / function evaluations
        double functionTol;                             // Function tolerance → stop when relative cost improvement < threshold
        double paramTol;                                // Parameter tolerance → stop when parameter updates are below threshold
        std::array<std::string_view, N> paramNames;     // Parameter names (for logging and diagnostics)
        double gradientTol{ 0.0 };                      // Gradient tolerance → stop when ∥Jᵀr∥ < threshold (stationarity)
        double lossScale{ 1.0 };                        // Loss function scale (e.g., δ for Huber or Cauchy)
        bool verbose{ false };                          // Logs the full Ceres calibration report
    };
}