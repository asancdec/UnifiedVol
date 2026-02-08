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
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Pricer.hpp"
#include "Optimization/Ceres/Optimizer.hpp"
#include "Optimization/Cost.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM
);

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM
);

template <std::floating_point T, std::size_t N>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Pricer<T, N>& pricer
);

namespace detail
{

template <std::floating_point T>
void validateInputs(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM
);

template <std::floating_point T, std::size_t N> struct PriceResidualJac;

std::array<double, 5> initGuess() noexcept;

std::array<double, 5> lowerBounds() noexcept;

std::array<double, 5> upperBounds() noexcept;

} // namespace detail
} // namespace uv::models::heston

#include "Models/Heston/Detail/Calibrate.inl"