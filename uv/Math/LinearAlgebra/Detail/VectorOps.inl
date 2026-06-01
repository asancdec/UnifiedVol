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

#include "Base/Macros/Require.hpp"
#include "Base/Types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace uv::math::linear_algebra
{
template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept
{
    evalInplace<T, N, F>(grid, std::forward<F>(f));
    return grid;
}

template <std::floating_point T, std::size_t N, typename F>
void evalInplace(std::array<T, N>& grid, F&& f) noexcept
{
    for (std::size_t i{0}; i < N; ++i)
        grid[i] = f(grid[i]);
}

template <std::floating_point T> T sum(std::span<const T> x) noexcept
{
    return std::accumulate(x.begin(), x.end(), T{});
}

template <std::floating_point T>
Vector<T> multiply(std::span<const T> v, const T x) noexcept
{
    std::size_t vSize{v.size()};

    Vector<T> result(vSize);

    for (std::size_t i = 0; i < vSize; ++i)
    {
        result[i] = v[i] * x;
    }

    return result;
}

template <std::floating_point T> Vector<T> reciprocal(std::span<const T> v) noexcept
{
    std::size_t vSize{v.size()};

    Vector<T> result(vSize);

    for (std::size_t i = 0; i < vSize; ++i)
    {
        result[i] = T(1.0) / v[i];
    }

    return result;
}

template <std::floating_point T>
Vector<T> hadamard(std::span<const T> a, std::span<const T> b)
{

    UV_REQUIRE_SAME_SIZE(a, b);

    std::size_t aSize{a.size()};

    Vector<T> c(aSize);

    for (std::size_t i = 0; i < aSize; ++i)
    {
        c[i] = a[i] * b[i];
    }

    return c;
}

template <typename T> Vector<T> makeSequence(std::size_t n, T start) noexcept
{
    Vector<T> v(n);
    std::iota(v.begin(), v.end(), start);
    return v;
}

template <typename T> T minValue(std::span<const T> x)
{
    UV_REQUIRE_NON_EMPTY(x);

    return *std::min_element(x.begin(), x.end());
}

template <typename T> T maxValue(std::span<const T> x)
{
    UV_REQUIRE_NON_EMPTY(x);

    return *std::max_element(x.begin(), x.end());
}

} // namespace uv::math::linear_algebra
