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

#include <cstddef>
#include <concepts>
#include <span>
#include <array>

namespace uv::math::pde
{
    /**
     * @brief Discrete initial density for Fokker–Planck evolution.
     *
     * Approximates a Dirac delta at x0 by assigning unit mass to the nearest node
     * (scaled by the local step size). Requires x0 within the grid domain.
     */
    template <std::floating_point T, std::size_t N>
    constexpr std::array<T, N> fokkerPlanckInit(T x0,
        std::span<const T, N> xGrid) noexcept;

    /**
     * @brief 1D Fokker–Planck solver (Chang–Cooper + Thomas).
     *
     * Marches `pdfGrid` forward `nT-1` implicit steps on a uniform grid, building a
     * Chang–Cooper tridiagonal system each step and solving it with the Thomas algorithm.
     *
     * Reference: Chang, J. S. & Cooper, G., "A Practical Difference Scheme for Fokker-Planck Equations",
     * Journal of Computational Physics 6, 1–16 (1970).
     */
    template<std::floating_point T, std::size_t nT, std::size_t nX>
    void fokkerPlanckSolve(std::span<T, nX> pdfGrid,
        std::span<T, nX - 1> B,
        std::span<T, nX - 1> C,
        T dt,
        T dx) noexcept;

    /**
     * @brief Solve a 1D Fokker–Planck system and log mass + runtime.
     *
     * Runs `fokkerPlanckSolve` for a log-space grid and prints a short diagnostic:
     * total mass (trapezoidal) and elapsed time, together with nT and nX.
     *
     * Reference: Chang, J. S. & Cooper, G., Journal of Computational Physics 6, 1–16 (1970).
     */
    template<std::floating_point T, std::size_t nT, std::size_t nX>
    void fokkerPlanckLog(std::span<T, nX> pdfGrid,
        std::span<T, nX - 1> B,
        std::span<T, nX - 1> C,
        T dt,
        T dx) noexcept;

    namespace detail
    {

        /**
         * @brief Compute a Chang–Cooper weight using a Taylor approximation.
         *
         * Uses omega = dx * B / C and a short Taylor series to avoid expensive exponentials.
         *
         * Reference: Chang, J. S. & Cooper, G., "A Practical Difference Scheme for Fokker-Planck Equations",
         * Journal of Computational Physics 6, 1–16 (1970).
         */
        template <std::floating_point T>
        constexpr T changCooperWeight(T B,
            T C,
            T dx) noexcept;

        /**
         * @brief Build Chang–Cooper tridiagonal coefficients (with reflecting boundaries) in one pass.
         *
         * Computes the tridiagonal diagonals `lower/middle/upper` for the Chang–Cooper
         * discretization of a 1D Fokker–Planck operator, and applies reflective (zero-flux)
         * boundary conditions at both ends of the grid.
         *
         * Weights are computed on-the-fly (rolling left/right) to avoid storing an
         * intermediate weights array.
         *
         * Reference: Chang, J. S. & Cooper, G., "A Practical Difference Scheme for Fokker-Planck Equations",
         * Journal of Computational Physics 6, 1–16 (1970).
         */
        template <std::floating_point T, std::size_t N>
        void changCooperDiagonals(std::span<T, N> upper,
            std::span<T, N> middle,
            std::span<T, N> lower,
            std::span<const T, N - 1> B,
            std::span<const T, N - 1> C,
            T dx,
            T invdx,
            T alpha) noexcept;

        /**
         * @brief Solve a tridiagonal linear system using the Thomas algorithm.
         *
         * Solves A * x = rhs in-place, where A is tridiagonal with:
         * - lower[i]  = sub-diagonal entry
         * - middle[i] = main diagonal entry b_i
         * - upper[i]  = super-diagonal entry c_i
         *
         * References:
         * - Tridiagonal matrix algorithm (Thomas algorithm), Wikipedia.
         *   https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
         */
        template <std::floating_point T, std::size_t N>
        void thomasSolve(std::span<T, N> x,
            std::span<const T, N> upper,
            std::span<const T, N> middle,
            std::span<const T, N> lower,
            std::span<T, N> scratch) noexcept;

    } // namespace detail

} // namespace uv::math::pde

#include "Functions.inl"