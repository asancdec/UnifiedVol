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

#include "Math/Integration/Functions.hpp"
#include "Utils/IO/Log.hpp"
#include "Utils/Aux/StopWatch.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <format>

namespace uv::math::pde
{
    template <std::floating_point T, std::size_t N>
    constexpr std::array<T, N> fokkerPlanckInit(T x0,
        std::span<const T, N> xGrid) noexcept
    {
        // ---------- Check size ----------

        static_assert
            (
                N >= 2, 
                "fokkerPlanckInit: grid must have at least 2 points"
             );

        // ---------- Binary search index ----------

        auto it = std::lower_bound(xGrid.begin(), xGrid.end(), x0);
        std::size_t idx{ static_cast<std::size_t>(it - xGrid.begin()) };

        // Clamp
        if (idx == N) idx = N - 1;

        // If there is an element to the left
        if (idx > 0)
        {
            // Choose the closest of idx and idx-1
            if (std::abs(xGrid[idx] - x0) > std::abs(xGrid[idx - 1] - x0))
            {
                --idx;
            }
        }

        // ---------- Dirac approximation ----------

        // Zero initialized
        std::array<T, N> p0{};
        T dX{};

        // Step size at the left 
        if (idx == 0)
        {
            dX = xGrid[1] - xGrid[0];
        }
        // Step size at the right
        else if (idx == N - 1)
        {
            dX = xGrid[N - 1] - xGrid[N - 2];
        }
        else
        {
            // Average step sizes
            const T dXLeft{ xGrid[idx] - xGrid[idx - 1] };
            const T dXRight{ xGrid[idx + 1] - xGrid[idx] };
            dX = T{ 0.5 } * (dXLeft + dXRight);
        }

        // Probability mass under p0 is 1
        // Therefore the following must hold approximately
        p0[idx] = T{ 1 } / dX;

        return p0;
    }

    template<std::floating_point T, std::size_t nT, std::size_t nX>
    void fokkerPlanckSolve(std::span<T, nX> pdfGrid,
        std::span<T, nX-1> B,
        std::span<T, nX-1> C,
        T dt,
        T dx) noexcept
    {
        // ---------- Check size ----------

        static_assert(nX >= 2, "fokkerPlanckSolve: nX must be >= 2");
        static_assert(nT >= 2, "fokkerPlanckSolve: nT must be >= 2");

        // ---------- Set ----------

        // Buffers
        std::array<T, nX> upper{};
        std::array<T, nX> middle{};
        std::array<T, nX> lower{};
        std::array<T, nX> scratch{};

        // Precompute
        T invdx{ T{ 1 } / dx };
        T alpha{ dt * invdx };

        for (std::size_t i = 0; i < nT - 1; ++i)
        {
            // ---------- Tridiagonal coefficients  ----------

            detail::changCooperDiagonals<T, nX>
                (
                    upper,
                    middle,
                    lower,
                    B,
                    C,
                    dx,
                    invdx,
                    alpha
                );

            // ---------- Solve ----------

            detail::thomasSolve<T, nX>
                (
                    pdfGrid,
                    upper,
                    middle,
                    lower,
                    scratch
                );
        }
    }

    template<std::floating_point T, std::size_t nT, std::size_t nX>
    void fokkerPlanckLog(std::span<T, nX> pdfGrid,
        std::span<T, nX - 1> B,
        std::span<T, nX - 1> C,
        T dt,
        T dx) noexcept
    {
        // Start timer
        utils::StopWatch timer_;
        timer_.StartStopWatch();

        // Solve
        math::pde::fokkerPlanckSolve<Real, nT, nX>
            (
                pdfGrid,
                B,
                C,
                dt,
                dx
            );

        // Stop timer
        timer_.StopStopWatch();

        // Log information
        UV_INFO
        (
            std::format
            (
                "[Fokker-Planck] mass={:.8e} ({:.4f} s nT={:L} nX={:L})",
                integration::trapezoidal(pdfGrid, dx),
                timer_.GetTime<std::ratio<1>>(),
                nT,
                nX
            )
        );
    }
} // namespace uv::math::pde

namespace uv::math::pde::detail
{
    template <std::floating_point T, std::size_t N>
    void changCooperDiagonals(std::span<T, N> upper,
        std::span<T, N> middle,
        std::span<T, N> lower,
        std::span<const T, N - 1> B,
        std::span<const T, N - 1> C,
        T dx,
        T invdx,
        T alpha) noexcept
    {

        // ---------- Dimension checks ----------

        static_assert
            (
                N >= 2,
                "changCooperDiagonals: grid must have at least 2 points"
                );

        // ---------- Calculate coefficients ----------

        // Precompute
        const T b0{ B[0] };
        const T c0{ C[0] };

        // Calculate the first weight
        T weight{ changCooperWeight(b0, c0, dx) };

        // Left boundary 
        upper[0] = -alpha * ((T{ 1 } - weight) * b0 + invdx * c0);
        middle[0] = T{ 1 } + alpha * (invdx * c0 - weight * b0);
        lower[0] = T{ 0 };

        for (std::size_t j = 1; j < N - 1; ++j)
        {
            // Calculate the next weight
            const T wRight{ changCooperWeight(B[j], C[j], dx) };

            // Precompute
            const T bRight{ B[j] };
            const T bLeft{ B[j - 1] };

            const T cRight{ C[j] };
            const T cLeft{ C[j - 1] };

            // Update coefficients
            upper[j] = -alpha *
                (
                    (T{ 1 } - wRight) * bRight
                    + invdx * cRight
                    );

            middle[j] = T{ 1 } + alpha *
                (
                    invdx * (cRight + cLeft)
                    + (T{ 1 } - weight) * bLeft
                    - wRight * bRight
                    );

            lower[j] = -alpha *
                (
                    invdx * cLeft
                    - weight * bLeft
                    );

            // Update the weight
            weight = wRight;
        }

        // Right boundary

        const std::size_t k{ N - 2 };
        const T bK{ B[k] };
        const T cK{ C[k] };

        upper[N - 1] = T{ 0 }; 
        middle[N - 1] = T{ 1 } + alpha * (invdx * cK + (T{ 1 } - weight) * bK);
        lower[N - 1] = -alpha * (invdx * cK - weight * bK);
    }

    template <std::floating_point T>
    constexpr T changCooperWeight(T B,
        T C,
        T dx) noexcept
    {
        // Precompute
        const T omega{ dx * B / C };
        const T omegaCubed{ omega * omega * omega };

        // Taylor series expansion
        return  T{ 0.5 } - omega / T{ 12 }
        + omegaCubed / T{ 720 }
        - omegaCubed * omega * omega / T{ 30240 };
    }

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