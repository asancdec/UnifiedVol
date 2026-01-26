// SPDX-License-Identifier: Apache-2.0
/*
 * File:        TanHSinH.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

namespace uv::math::integration
{
/**
 * @brief Fixed Tanh–Sinh quadrature rule for integrals on (0, +∞).
 *
 * Implements a precomputed (compile-time sized) Tanh–Sinh grid and provides:
 * - scalar integration on (0, +∞)
 * - multi-output integration (vectorised integrand) on (0, +∞)
 * - grid printing for debugging/inspection
 *
 * The rule uses a change of variables mapping (0, +∞) → (-1, 1) with
 * Tanh–Sinh node generation and stores transformed inputs and weights.
 *
 * @tparam N Number of quadrature nodes (must be even).
 *
 * @note This class is header-only (see TanHSinH.inl).
 */
template <std::floating_point T, std::size_t N> class TanHSinH
{
  private:
    //--------------------------------------------------------------------------
    // Internal struct
    //--------------------------------------------------------------------------
    /**
     * @brief Precomputed node data for the Tanh–Sinh rule.
     *
     * Stores the raw node/weight (x,w) plus precomputed transformed inputs and
     * scaling factors for the right-hand side (RHS) and left-hand side (LHS)
     * contributions used in the (0, +∞) mapping.
     */
    struct Node
    {
        T w;           // Weight value
        T y;           // y_n term
        T x;           // Abscissa value
        T factorRight; // Scaling factor (RHS)
        T inputRight;  // Transformed input (RHS)
        T factorLeft;  // Scaling factor (LHS)
        T inputLeft;   // Transformed input (LHS)
    };

    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------

    const T h_;                 // Step size
    std::array<Node, N> nodes_; // Precomputed node storage

    //--------------------------------------------------------------------------
    // Math
    //--------------------------------------------------------------------------

    /**
     * @brief Generate a single precomputed node for a given n*h.
     *
     * @param nh Value of n*h used to build the node.
     * @return Fully populated Node (weights + transformed inputs).
     */
    Node generateNode(T nh) const noexcept;

  public:
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    /**
     * @brief Construct and precompute the Tanh–Sinh nodes.
     *
     * Computes an “optimal” step size (heuristic) and fills @ref nodes_.
     *
     * @note Compile-time checks typically enforce:
     * - N > 0
     * - N is even (for unroll-by-2 loops)
     */
    TanHSinH();

    //--------------------------------------------------------------------------
    // Math
    //--------------------------------------------------------------------------

    /**
     * @brief Numerically integrate a scalar-valued function on (0, +∞).
     *
     * @tparam F Callable type with signature: T f(T x).
     * @param f Integrand callable.
     * @return Numerical approximation of ∫₀^∞ f(x) dx.
     *
     * @note Uses early-exit checks when additional terms are negligible.
     */
    template <typename F> T integrateZeroToInf(F&& f) const noexcept;

    /**
     * @brief Numerically integrate a multi-output function on (0, +∞).
     *
     * Computes M integrals in one pass, useful when the integrand returns
     * multiple components (e.g. vector of payoffs).
     *
     * @tparam M Number of components returned by the integrand.
     * @tparam F Callable type with signature: std::array<T, M> f(T x).
     * @param f Integrand callable.
     * @return Array of M integrals, component-wise.
     *
     * @note Uses component-wise early-exit checks per accumulated sum.
     */
    template <std::size_t M, typename F>
    std::array<T, M> integrateZeroToInfMulti(F&& f) const noexcept;

    //--------------------------------------------------------------------------
    // Utilities
    //--------------------------------------------------------------------------

    /**
     * @brief Print the precomputed grid (nodes and weights) to the logger.
     *
     * Intended for debugging and verification of the quadrature grid.
     */
    void printGrid() const noexcept;
};
} // namespace uv::math::integration

#include "TanHSinH.inl"
