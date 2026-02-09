// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

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

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::models::heston::calibrate
{
template <
    std::floating_point T,
    std::size_t N,
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
    const core::Matrix<T>& callM,
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
    const core::Matrix<double>& callM,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
);

template <std::floating_point T>
void validateInputs(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM
);

std::array<double, 5> initGuess() noexcept;

std::array<double, 5> lowerBounds() noexcept;

std::array<double, 5> upperBounds() noexcept;

} // namespace detail
} // namespace uv::models::heston::calibrate

#include "Models/Heston/Calibrate/Detail/Calibrate.inl"