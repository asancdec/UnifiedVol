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

#include "Core/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/Aux/StopWatch.hpp"
#include "Utils/IO/Log.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <iterator>
#include <string>
#include <utility>

namespace uv::math::pde
{
template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> andreasenHugeInit(const std::array<T, N>& xGrid, F&& payoff)
{

    auto payoffLogSpace = [payoff](T x) -> T
    {
        return payoff(1.0, std::exp(x));
    };

    return core::eval(xGrid, payoffLogSpace);
}

template <std::floating_point T, std::size_t NT, std::size_t NX>
void andreasenHugeSolve(
    std::span<T, NX> c,
    std::span<T, NX - 2> cInner,
    AHCache<T, NX>& aHCache
) noexcept
{
    static_assert(NX >= 4, "Needs at least 4 grid points for boundaries.");

    constexpr std::size_t innerNX{NX - 2};

    const T aL{aHCache.aL};
    const T bL{aHCache.bL};
    const T aR{aHCache.aR};
    const T bR{aHCache.bR};

    const T zLower{aHCache.zLower};
    const T zMiddle{aHCache.zMiddle};
    const T zUpper{aHCache.zUpper};

    std::array<T, NX - 2>& scratch{aHCache.scratch};
    std::array<T, NX - 2>& lower{aHCache.lower};
    std::array<T, NX - 2>& middle{aHCache.middle};
    std::array<T, NX - 2>& upper{aHCache.upper};

    const std::array<T, NX - 2>& localVar{aHCache.localVar};

    const T& lowerFirst{lower.front()};
    T& middleFirst{middle.front()};
    T& upperFirst{upper.front()};

    T& lowerLast{lower.back()};
    T& middleLast{middle.back()};
    const T& upperLast{upper.back()};

    for (std::size_t i{0}; i < NT; ++i)
    {

        for (std::size_t j{0}; j < innerNX; ++j)
        {
            const T localVarCurr{localVar[j]};

            lower[j] = zLower * localVarCurr;
            middle[j] = zMiddle * localVarCurr + 1.0;
            upper[j] = zUpper * localVarCurr;
        }

        middleFirst += lowerFirst * aL;
        upperFirst += lowerFirst * bL;

        lowerLast += upperLast * bR;
        middleLast += upperLast * aR;

        detail::thomasSolve<T, NX - 2>(cInner, upper, middle, lower, scratch);
    }

    c.front() = aL * c[1] + bL * c[2];
    c.back() = aR * c[NX - 2] + bR * c[NX - 3];
}
} // namespace uv::math::pde

namespace uv::math::pde::detail
{
template <std::floating_point T, std::size_t N>
void thomasSolve(
    std::span<T, N> x,
    std::span<const T, N> upper,
    std::span<const T, N> middle,
    std::span<const T, N> lower,
    std::span<T, N> scratch
) noexcept
{

    static_assert(N >= 2, "thomasSolve: N must be >= 2");

    scratch[0] = upper[0] / middle[0];
    x[0] = x[0] / middle[0];

    for (std::size_t i = 1; i < N; i++)
    {

        const T lowerI{lower[i]};
        const T invDenom{1.0 / (middle[i] - lowerI * scratch[i - 1])};

        if (i < N - 1)
        {
            scratch[i] = upper[i] * invDenom;
        }
        x[i] = (x[i] - lowerI * x[i - 1]) * invDenom;
    }

    for (std::size_t i = N - 1; i-- > 0;)
    {
        x[i] -= scratch[i] * x[i + 1];
    }
}
} // namespace uv::math::pde::detail