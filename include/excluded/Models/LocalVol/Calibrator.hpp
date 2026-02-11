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

#include "Core/Matrix/Matrix.hpp"
#include "Core/Types.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Optimization/Ceres/Optimizer.hpp"
#include "Math/Optimization/Functions.hpp"
#include "Models/LocalVol/Pricer.hpp"
#include "Models/LocalVol/Surface.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::models::localvol::calibrator
{

namespace opt = uv::math::opt;

template <
    std::floating_point T,
    std::size_t NT,
    std::size_t NX,
    class Interpolator,
    typename Policy>
Surface<T> calibrate(
    const core::VolSurface<T>& volSurface,
    Pricer<T, NT, NX, Interpolator>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
);

template <
    std::floating_point T,
    std::size_t NT,
    std::size_t NX,
    class Interpolator,
    typename Policy>
Surface<T> calibrate(
    const core::Matrix<T>& callPrice,
    const Vector<T>& maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totVar,
    Pricer<T, NT, NX, Interpolator>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
);

namespace detail
{
template <std::floating_point T> void validate(
    const core::Matrix<T>& callPrice,
    const Vector<T>& maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totVar
);

template <std::floating_point T>
Vector<double> coldGuess(T tenor, std::span<const T> totVar) noexcept;

template <std::floating_point T>
void warmGuess(std::span<const T> params, std::span<T> localVar) noexcept;

} // namespace detail
} // namespace uv::models::localvol::calibrator

#include "Calibrator.inl"