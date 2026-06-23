// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Optimization/Ceres/Optimizer.hpp"

namespace uv::models::heston::calibrate::detail
{

template <typename Policy> void setGuessBounds(opt::ceres::Optimizer<Policy>& optimizer);

} // namespace uv::models::heston::calibrate::detail

#include "Models/Heston/Calibrate/Detail/Initialize.inl"
