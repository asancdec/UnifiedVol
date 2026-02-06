// SPDX-License-Identifier: Apache-2.0
/*
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

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace uv::opt
{
void clampBounds(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate = true
);

void clampLowerBounds(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    bool doValidate = true
);

void clampUpperBounds(
    std::span<double> initGuess,
    std::span<const double> upperBounds,
    bool doValidate = true
);

void warnBoundsHit(
    std::span<const double> x,
    const std::optional<std::vector<double>>& lowerBounds,
    const std::optional<std::vector<double>>& upperBounds,
    bool doValidate = true
);

void warnBoundsHit(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate = true
);

void logResults(
    std::span<const double> x,
    std::span<const std::string_view> paramNames,
    double sse,
    unsigned iterCount,
    double elapsedMs,
    bool isSuccess
);

void validateBounds(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
);

void validateBoundsSpec(
    std::size_t n,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
);

void validateLowerBounds(std::span<const double> x, std::span<const double> lowerBounds);

void validateLowerBoundsSpec(std::size_t n, std::span<const double> lowerBounds);

void validateUpperBounds(std::span<const double> x, std::span<const double> upperBounds);

void validateUpperBoundsSpec(std::size_t n, std::span<const double> upperBounds);
} // namespace uv::opt
