// SPDX-License-Identifier: Apache-2.0
/*
 * File:        AHCache.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-20
 *
 * Description:
 *   Cache and reusable buffers for the Andreasen–Huge (AH) forward PDE solver
 *   on a uniform x-grid (log-moneyness / log-forward-moneyness).
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

namespace uv::math::pde
{
//--------------------------------------------------------------------------
/**
 * @brief Cached coefficients and work buffers for the Andreasen–Huge PDE step.
 *
 * @details
 * Internal cache used by the Andreasen–Huge solver to avoid recomputing
 * grid-dependent constants and to reuse tridiagonal work buffers.
 *
 * The normalized call price C(x,t) is advanced on a uniform x-grid using a
 * Crank–Nicolson scheme, resulting in a tridiagonal system solved in-place
 * via a Thomas algorithm.
 *
 * Simple linear boundary closures are applied at the first and last grid
 * points:
 *
 *   C_0     = aL*C_1     + bL*C_2
 *   C_{N-1} = aR*C_{N-2} + bR*C_{N-3}
 *
 * Grid-dependent terms are computed once in the pricer constructor, while
 * time-step-dependent coefficients are updated per maturity.
 */
template <std::floating_point T, std::size_t N> struct AHCache
{
    // ---------- Grid-dependent constants (set once) ----------

    T invDXSquared{};         ///< 1 / dx^2, where dx is the uniform x-grid spacing.
    T invDXSquaredTimesTwo{}; ///< 2 / dx^2 (precomputed to avoid repeated
                              ///< multiply).

    T lowerFactor{}; ///< (1/dx^2 + 1/(2dx)) factor used in lower diagonal
                     ///< assembly.
    T upperFactor{}; ///< (1/dx^2 - 1/(2dx)) factor used in upper diagonal
                     ///< assembly.

    // ---------- Boundary closure coefficients (set once) ----------

    T aL{}; ///< Left boundary: C_0 = aL*C_1 + bL*C_2.
    T bL{}; ///< Left boundary: C_0 = aL*C_1 + bL*C_2.
    T aR{}; ///< Right boundary: C_{N-1} = aR*C_{N-2} + bR*C_{N-3}.
    T bR{}; ///< Right boundary: C_{N-1} = aR*C_{N-2} + bR*C_{N-3}.

    // ---------- Time-step coefficients (update per dt) ----------

    T zLower{};  ///< zLower  = -(dt/2) * lowerFactor.
    T zMiddle{}; ///< zMiddle =  (dt/2) * (2/dx^2).
    T zUpper{};  ///< zUpper  = -(dt/2) * upperFactor.

    // ---------- Work buffers (reused every time step) ----------

    std::array<T, N - 2>
        scratch{}; ///< Thomas scratch buffer (e.g. modified upper or temp).
    std::array<T, N - 2> lower{};  ///< Tridiagonal lower diagonal (interior nodes only).
    std::array<T, N - 2> middle{}; ///< Tridiagonal main diagonal (interior nodes only).
    std::array<T, N - 2> upper{};  ///< Tridiagonal upper diagonal (interior nodes only).

    std::array<T, N - 2> localVar{}; ///< Local variance grid
};
} // namespace uv::math::pde
