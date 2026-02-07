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

#include "Models/SVI/Calibrate/Detail/SliceData.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <array>
#include <concepts>

namespace uv::models::svi::detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo>
void setGuessBounds(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const Params<T>* prevParams,
    const SliceData& sliceData
) noexcept;

std::array<double, 4> coldGuess() noexcept;

std::array<double, 4> warmGuess(const Params<double>& params) noexcept;

std::array<double, 4> lowerBounds(double logKFMin) noexcept;

std::array<double, 4> upperBounds(double logKFMax) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Initialize.inl"