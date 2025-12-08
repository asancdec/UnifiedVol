// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include <algorithm>
#include <string> 
#include <array>
#include <cstddef>
#include <string_view>
#include <span>

namespace uv::math::opt
{
    // Clamp initial guess within upper and lower bounds
    template <std::size_t N>
    void clamp(std::array<Real, N>& initGuess,
        const std::array<Real, N>& lowerBounds,
        const std::array<Real, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Warn if upper or lower bounds are touched
    template <std::size_t N>
    void warnBoundsHit(std::span<Real> x,
        const std::array<Real, N>& lowerBounds,
        const std::array<Real, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Log calibration results 
    template <std::size_t N>
    void logResults(std::span<Real> x,
        const std::array<std::string_view, N>& paramNames,
        Real sse,
        unsigned iterCount,
        Real elapsedMs,
        bool isSuccess) noexcept;
}

#include "Functions.inl"