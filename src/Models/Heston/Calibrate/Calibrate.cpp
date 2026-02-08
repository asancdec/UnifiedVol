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

#include "Models/Heston/Calibrate/Calibrate.hpp"

namespace uv::models::heston::detail
{
std::array<double, 5> initGuess() noexcept
{
    return {2.5, 0.09, 0.60, -0.75, 0.09};
}

std::array<double, 5> lowerBounds() noexcept
{
    return {0.001, 0.001, 0.001, -0.999, 0.001};
}

std::array<double, 5> upperBounds() noexcept
{
    return {10.0, 0.5, 10.0, 0.999, 0.5};
}
} // namespace uv::models::heston::detail
