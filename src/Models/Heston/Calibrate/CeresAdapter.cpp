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

#include "Models/Heston/Calibrate/CeresAdapter.hpp"

namespace uv::models::heston::calibrate::detail
{
opt::ceres::Config makeCeresConfig(const Config& config) noexcept
{

    return opt::ceres::Config{
        .maxEval = config.maxEval,
        .functionTol = config.tolerance,
        .paramTol = config.tolerance,
        .gradientTol = config.tolerance,
        .paramNames = paramNames,
        .verbosity = config.verbosity,
        .numThreads = config.numThreads
    };
}

opt::ceres::Optimizer<HestonPolicy> makeOptimizer(const Config& config) noexcept
{
    return opt::ceres::Optimizer<HestonPolicy>{makeCeresConfig(config)};
}

} // namespace uv::models::heston::calibrate::detail