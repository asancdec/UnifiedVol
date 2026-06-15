// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Calibrate/CeresAdapter.hpp"

namespace uv::models::heston::calibrate::detail
{
opt::ceres::Config makeCeresConfig(const Config& config)
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

opt::ceres::Optimizer<HestonPolicy> makeOptimizer(const Config& config)
{
    return opt::ceres::Optimizer<HestonPolicy>{makeCeresConfig(config)};
}

} // namespace uv::models::heston::calibrate::detail
