// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Price/Pricer.hpp"
#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Optimizer.hpp"
#include "Optimization/Cost.hpp"

#include <concepts>
#include <cstddef>
#include <span>

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

namespace detail
{
template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
);

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& vol,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
);

template <
    std::floating_point CalcT,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<double> calibrateDouble(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& vol,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
);

} // namespace detail
} // namespace uv::models::heston::calibrate

#include "Models/Heston/Calibrate/Detail/Calibrate.inl"