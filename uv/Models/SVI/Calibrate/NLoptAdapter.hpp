// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "Models/SVI/Calibrate/Config.hpp"
#include "Optimization/NLopt/Config.hpp"

namespace uv::models::svi::detail
{
inline opt::nlopt::Config<4> makeNLoptConfig(const Config& config) noexcept
{

    return opt::nlopt::Config<4>{
        .tol = config.objectiveTol,
        .ftolRel = config.objectiveTol,
        .maxEval = config.maxEval,
        .verbose = config.verbose,
        .paramNames = paramNames
    };
}

} // namespace uv::models::svi::detail