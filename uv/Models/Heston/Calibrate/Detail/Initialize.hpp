// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Optimization/Ceres/Optimizer.hpp"

#include <array>

namespace uv::models::heston::calibrate::detail
{

template <typename Policy> void setGuessBounds(opt::ceres::Optimizer<Policy>& optimizer);

constexpr std::array<double, 5> initGuess() noexcept
{
    return {2.5, 0.09, 1.0, -0.75, 0.20};
}

constexpr std::array<double, 5> lowerBounds() noexcept
{
    return {0.001, 0.001, 0.001, -0.999, 0.001};
}

constexpr std::array<double, 5> upperBounds() noexcept
{
    return {15.0, 0.5, 10.0, 0.999, 0.5};
}

} // namespace uv::models::heston::calibrate::detail

#include "Models/Heston/Calibrate/Detail/Initialize.inl"