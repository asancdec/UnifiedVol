// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Config.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Ceres optimizer parameters configuration
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

namespace uv::math::opt::ceres
{
//--------------------------------------------------------------------------
/**
 * @brief Configuration for Ceres-based optimisation.
 *
 * @details
 * Groups numerical tolerances, iteration limits, robust loss scaling,
 * and diagnostic options used by the Ceres optimiser wrapper.
 *
 * Solver strategy (trust-region method, linear solver, loss type) is
 * selected at compile time via the optimisation Policy, while this
 * structure provides runtime tuning parameters.
 *
 * @tparam N Number of optimisation parameters.
 */
template <std::size_t N> struct Config
{
    // ---------- Termination ----------

    unsigned maxEval;   ///< Maximum number of function evaluations.
    double functionTol; ///< Stop when relative cost reduction falls below this
                        ///< threshold.
    double paramTol;    ///< Stop when parameter updates fall below this threshold.
    double gradientTol; ///< Stop when norm of J^T r falls below this threshold.

    // ---------- Robust loss  ----------

    double lossScale; ///< Robust loss scale (e.g. delta for Huber), in residual
                      ///< units.

    // ---------- Logging ----------

    std::array<std::string_view, N>
        paramNames;      ///< Parameter names for logging and diagnostics.
    bool verbose{false}; ///< Enable full Ceres solver report output.
};
} // namespace uv::math::opt::ceres