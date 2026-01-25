// SPDX-License-Identifier: Apache-2.0
/*
 * File:        VarianceView.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-01-23
 *
 * Description:
 *   Lightweight non-owning view over a single local-variance slice.
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

#include <span>
#include <concepts>

namespace uv::models::localvol
{
    /**
     * @brief Non-owning view of a local-variance row at a fixed maturity.
     *
     * Provides read-only access to the interpolation grid (log-forward moneyness),
     * local variance values, and precomputed derivatives used by interpolation
     * and PDE solvers.
     *
     * All spans are required to have the same length and remain valid for the
     * lifetime of the view.
     */
    template <std::floating_point T>
    struct VarianceView
    {
        std::span<const T> xs;    ///< log(F/K) grid nodes.
        std::span<const T> ys;    ///< Local variance values at grid nodes.
        std::span<const T> dydx;  ///< Derivatives d(var)/d(log(F/K)) at nodes.
    };
} // namespace uv::models::localvol