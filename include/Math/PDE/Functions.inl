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

#include "Utils/IO/Log.hpp"
#include "Utils/Aux/StopWatch.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Core/Functions.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <format>

namespace uv::math::pde
{   
    template 
    <
        std::floating_point T, 
        std::size_t N,
        typename F
    >
        std::array<T, N> andreasenHugeInit(const std::array<T, N>& xGrid, F&& payoff)
    {
        // ---------- Bind to log-space ----------

        auto payoffLogSpace = [payoff](T x) -> T
        {
            return payoff(1.0, std::exp(x));
        };

        // ---------- Evaluate payoff function ----------

        return core::eval
        (
            xGrid, 
            payoffLogSpace
        );
    }

    template
        <
        std::floating_point T,
        std::size_t NT,
        std::size_t NX
        >
        void andreasenHugeSolve(std::array<T, NX>& c,
            const std::array<T, NX>& cInit,
            const std::array<T, NX>& localVar,
            models::localvol::AHCache<T, NX>& aHCache)
    {

        // Write the initial values into the buffer
        c = cInit;

        //for (std::size_t i{ 0 }; i < NT; ++i)
        //{
        //    lower = localVar;

        //}
    }



} // namespace uv::math::pde

namespace uv::math::pde::detail
{
    template <std::floating_point T, std::size_t N>
    void thomasSolve(std::span<T, N> x,
        std::span<const T, N> upper,
        std::span<const T, N> middle,
        std::span<const T, N> lower,
        std::span<T, N> scratch) noexcept
    {
        // ---------- Dimension checks ----------

        static_assert
        (
            N >= 2,
            "thomasSolve: N must be >= 2"
        );

        // ---------- Solve ----------

        scratch[0] = upper[0] / middle[0];
        x[0] = x[0] / middle[0];

        for (std::size_t i = 1; i < N; i++)
        {
            // Precompute
            const T lowerI{ lower[i] };
            const T invDenom{ T{1} / (middle[i] - lowerI * scratch[i - 1]) };

            if (i < N - 1)
            {
                scratch[i] = upper[i] * invDenom;
            }
            x[i] = (x[i] - lowerI * x[i - 1]) * invDenom;
        }

        for (std::size_t i = N - 1; i-- > 0; )
        {
            x[i] -= scratch[i] * x[i + 1];
        }
    }
} // namespace uv::math::pde::detail