// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
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

#include "Utils/Aux/Errors.hpp"  
#include "Core/Functions.hpp"

#include <string> 
#include <vector>
#include <numeric>  
#include <algorithm>

namespace uv::core
{       
    template <std::floating_point T, std::size_t N>
    constexpr std::array<T, N> generateGrid(T bound1,
        T bound2) noexcept
    {
        // Check dimensions
        static_assert(N >= 2, "generateGrid: grid must have at least 2 points");

        std::array<T, N> grid{};
        grid[0] = bound1;

        // Calculate step size
        const T dx{ (bound2 - bound1) / static_cast<T>(N - 1) };

        for (std::size_t i = 1; i < N; ++i)
        {
            // Set grid elements
            grid[i] = bound1 + dx * static_cast<T>(i);
        }

        return grid;
    }

    template <std::floating_point T, std::size_t N, typename F>
    std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept
    {
        evalInplace<T, N, F>(grid, f);
        return grid;
    }

    template <std::floating_point T, std::size_t N, typename F>
    void evalInplace(std::array<T, N>& grid, F&& f) noexcept
    {
        for (std::size_t i{ 0 }; i < N; ++i)
            grid[i] = f(grid[i]);
    }

    template <std::floating_point T>
    T sum(std::span<const T> x) noexcept
    {
        return std::accumulate(x.begin(), x.end(), T{});
    }

    template <std::floating_point T>
    Vector<T> multiply(std::span<const T> v,
        const T x) noexcept
    {
        std::size_t vSize{ v.size() };

        Vector<T> result(vSize);

        for (std::size_t i = 0; i < vSize; ++i)
        {
            result[i] = v[i] * x;
        }

        return result;
    }

    template <std::floating_point T>
    Vector<T> reciprocal(std::span<const T> v) noexcept
    {
        std::size_t vSize{ v.size() };

        Vector<T> result(vSize);

        for (std::size_t i = 0; i < vSize; ++i)
        {
            result[i] = T(1.0) / v[i];
        }

        return result;
    }

    template <std::floating_point T>
    Vector<T> hadamard(std::span<const T> a,
        std::span<const T> b)
    {
        // ---------- Check dimensions ----------

        std::size_t aSize{ a.size() };
        std::size_t bSize{ b.size() };

        UV_REQUIRE(
            aSize == bSize,
            ErrorCode::InvalidArgument,
            "hadamard: vectors must have same size (got " +
            std::to_string(aSize) + " and " + std::to_string(bSize) + ")"
        );

        // ---------- Element-wise multiplication ----------

        Vector<T> c(aSize);

        for (std::size_t i = 0; i < aSize; ++i)
        {
            c[i] = a[i] * b[i];
        }

        return c;
    }

    template <typename T>
    Vector<T> makeSequence(std::size_t n,
        T start) noexcept
    {
        Vector<T> v(n);
        std::iota(v.begin(), v.end(), start);
        return v;
    }

    template<typename T>
    T minValue(std::span<const T> x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "minValue: input vector is empty"
        );

        return *std::min_element(x.begin(), x.end());
    }

    template<typename T>
    T maxValue(std::span<const T> x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "maxValue: input vector is empty"
        );

        return *std::max_element(x.begin(), x.end());
    }


    template <typename To, typename From>
    Vector<To> convertVector(const Vector<From>& x) noexcept
    {
        Vector<To> out;
        out.reserve(x.size());

        for (const From& v : x)
            out.push_back(static_cast<To>(v));

        return out;
    }

} // namespace uv::core