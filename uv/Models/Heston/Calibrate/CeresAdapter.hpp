// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/Heston/Calibrate/Config.hpp"
#include "Optimization/Ceres/Optimizer.hpp"

namespace uv::models::heston::calibrate::detail
{
opt::ceres::Config makeCeresConfig(const Config& config) noexcept;

opt::ceres::Optimizer<HestonPolicy> makeOptimizer(const Config& config) noexcept;

} // namespace uv::models::heston::calibrate::detail