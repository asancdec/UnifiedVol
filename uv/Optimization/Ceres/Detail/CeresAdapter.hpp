// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Optimization/Ceres/Config.hpp"

#include <ceres/ceres.h>
#include <memory>

namespace uv::opt::ceres::detail
{
::ceres::TrustRegionStrategyType toCeres(TrustRegionStrategy a);

::ceres::LinearSolverType toCeres(LinearSolver a);

std::unique_ptr<::ceres::LossFunction> makeLoss(Loss a, double lossParam);

} // namespace uv::opt::ceres::detail