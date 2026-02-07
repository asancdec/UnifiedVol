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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "Base/Alias.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <chrono>
#include <concepts>
#include <iostream>
#include <span>

namespace uv::models::svi
{
template <std::floating_point T, opt::nlopt::Algorithm Algo>
Vector<Params<T>> calibrate(
    const core::VolSurface<T>& volSurface,
    const opt::nlopt::Optimizer<4, Algo>& prototype
);

template <std::floating_point T, opt::nlopt::Algorithm Algo>
Vector<Params<T>> calibrate(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype
);

namespace detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo>
Params<T> calibrateSlice(
    T t,
    std::span<const double> logKF,
    std::span<const double> totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    const Params<T>* prevParams
);

template <std::floating_point T>
void validateInputs(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance
);

} // namespace detail
} // namespace uv::models::svi

#include "Models/SVI/Calibrate/Detail/Calibrate.inl"