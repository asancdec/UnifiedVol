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

#include "Math/PDE/AHCache.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::pde
{

template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> andreasenHugeInit(const std::array<T, N>& xGrid, F&& payoff);

template <std::floating_point T, std::size_t NT, std::size_t NX>
[[gnu::hot]] void andreasenHugeSolve(
    std::span<T, NX> c,
    std::span<T, NX - 2> cInner,
    AHCache<T, NX>& aHCache
) noexcept;

namespace detail
{

template <std::floating_point T, std::size_t N>
[[gnu::hot]] void thomasSolve(
    std::span<T, N> x,
    std::span<const T, N> upper,
    std::span<const T, N> middle,
    std::span<const T, N> lower,
    std::span<T, N> scratch
) noexcept;

}

} // namespace uv::math::pde

#include "Functions.inl"