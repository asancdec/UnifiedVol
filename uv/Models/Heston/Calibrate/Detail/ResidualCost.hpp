// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/Heston/Calibrate/Detail/MaturitySlice.hpp"
#include "Models/Heston/Price/Pricer.hpp"
#include "Optimization/Ceres/Config.hpp"

#include <concepts>
#include <cstddef>
#include <memory>

#include <ceres/cost_function.h>

namespace uv::models::heston::calibrate::detail
{

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction>
makeSliceCost(const MaturitySlice& slice, const price::Pricer<T, N>& pricer);
} // namespace uv::models::heston::calibrate::detail

#include "Models/Heston/Calibrate/Detail/ResidualCost.inl"
