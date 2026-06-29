// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Optimization/NLopt/Config.hpp"

namespace uv::models::svi
{

struct Config
{
    double objectiveTol{1e-12};
    unsigned int maxEval{10000};
    bool verbose{true};
    bool printParams{false};
};

namespace detail
{

constexpr opt::nlopt::Config<4> makeNLoptConfig(const Config& config) noexcept
{
    return opt::nlopt::Config<4>{
        .tol = config.objectiveTol,
        .ftolRel = config.objectiveTol,
        .maxEval = config.maxEval,
        .verbose = config.verbose,
        .paramNames = {"b", "rho", "m", "sigma"}
    };
}

} // namespace detail

} // namespace uv::models::svi
