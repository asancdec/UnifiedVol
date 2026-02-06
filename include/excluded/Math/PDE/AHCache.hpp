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

#include <array>
#include <concepts>
#include <cstddef>

namespace uv::math::pde
{

template <std::floating_point T, std::size_t N> struct AHCache
{

    T invDXSquared;
    T invDXSquaredTimesTwo;

    T lowerFactor;

    T upperFactor;

    T aL;
    T bL;
    T aR;
    T bR;

    T zLower;
    T zMiddle;
    T zUpper;

    std::array<T, N - 2> scratch;
    std::array<T, N - 2> lower;
    std::array<T, N - 2> middle;
    std::array<T, N - 2> upper;

    std::array<T, N - 2> localVar;
};
} // namespace uv::math::pde
