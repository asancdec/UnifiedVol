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

#include <algorithm>
#include <string> 
#include <array>
#include <cstddef>
#include <string_view>
#include <span>

namespace uv::math::opt
{
    /**
     * @brief Clamp an initial parameter guess inside its admissible bounds.
     *
     * Each parameter is clamped component-wise:
     *   initGuess[i] ← clamp(initGuess[i], lowerBounds[i], upperBounds[i])
     *
     * A warning is logged whenever a parameter is modified.
     *
     * Typical use case:
     *   - Sanitising user-provided initial guesses before starting an optimiser.
     *
     * @tparam N Number of parameters.
     * @param[in,out] initGuess Initial parameter vector (modified in-place).
     * @param[in] lowerBounds Component-wise lower bounds.
     * @param[in] upperBounds Component-wise upper bounds.
     * @param[in] paramNames Human-readable parameter names (for logging).
     *
     * @note This function never throws and is safe to call inside optimisation setup.
     */
    template <std::size_t N>
    void clamp(std::array<double, N>& initGuess,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    /**
     * @brief Emit warnings if any parameter lies (numerically) on its bounds.
     *
     * A parameter is considered to have hit a bound if it is within
     * a small absolute/relative tolerance of the bound value.
     *
     * Typical use case:
     *   - Post-optimisation diagnostics to detect constrained solutions.
     *
     * @tparam N Number of parameters.
     * @param[in] x Optimised parameter vector.
     * @param[in] lowerBounds Component-wise lower bounds.
     * @param[in] upperBounds Component-wise upper bounds.
     * @param[in] paramNames Human-readable parameter names (for logging).
     *
     * @note This function does not modify the parameter vector.
     * @note Warnings are informational only and do not imply optimisation failure.
     */
    template <std::size_t N>
    void warnBoundsHit(std::span<double> x,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    /**
     * @brief Log optimisation results in a standardised calibration format.
     *
     * The output includes:
     *   - parameter values
     *   - sum of squared errors (SSE)
     *   - elapsed time (ms)
     *   - iteration count
     *   - success flag
     *
     * Typical output:
     *   [Calib] kappa=1.2345  theta=0.0456  ...  SSE=1.23e-06 (2.31 ms, 54 it, SUCCESS)
     *
     * @tparam N Number of parameters.
     * @param[in] x Optimised parameter vector.
     * @param[in] paramNames Human-readable parameter names.
     * @param[in] sse Final sum of squared errors.
     * @param[in] iterCount Number of optimiser iterations.
     * @param[in] elapsedMs Elapsed wall-clock time in milliseconds.
     * @param[in] isSuccess Whether the optimiser reported success.
     */
    template <std::size_t N>
    void logResults(std::span<double> x,
        const std::array<std::string_view, N>& paramNames,
        double sse,
        unsigned iterCount,
        double elapsedMs,
        bool isSuccess) noexcept;
}

#include "Functions.inl"