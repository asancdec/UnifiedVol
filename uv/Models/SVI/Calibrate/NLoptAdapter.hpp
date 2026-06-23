// SPDX-License-Identifier: Apache-2.0

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
        .paramNames = {"b", "rho", "m", "sigma"}
    };
}

} // namespace uv::models::svi::detail
