// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Price/Pricer.hpp"
#include "Optimization/Ceres/Config.hpp"

#include <concepts>
#include <cstddef>

namespace uv::models::heston::calibrate
{
template <
    std::floating_point T,
    std::size_t N = defaultNodes,
    opt::ceres::GradientMode Mode = HestonGradient,
    typename Policy = HestonPolicy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config = {}
);

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode = HestonGradient,
    typename Policy = HestonPolicy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config,
    price::Pricer<T, N>& pricer
);
} // namespace uv::models::heston::calibrate

#include "Models/Heston/Calibrate/Detail/Calibrate.inl"
