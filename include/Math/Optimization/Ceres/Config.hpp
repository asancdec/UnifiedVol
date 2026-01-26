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

#include <vector>
#include <array>
#include <cstddef>
#include <string_view>

namespace uv::math::opt::ceres
{

/**
 * @brief Configuration parameters for a Ceres-based optimisation.
 */
struct Config
{
    // ---------- Termination ----------

    unsigned maxEval;   ///< Maximum number of function evaluations.
    double functionTol; ///< Relative cost reduction tolerance.
    double paramTol;    ///< Parameter update tolerance.
    double gradientTol; ///< Gradient norm tolerance (|J^2 * r|).

    // ---------- Robust loss ----------

    double lossScale; ///< Robust loss scale.

    // ---------- Logging ----------

    /// Parameter names used for logging
    std::vector<std::string_view> paramNames;

    bool verbose{false}; ///< Enable detailed solver output.
};

} // namespace uv::math::opt::ceres