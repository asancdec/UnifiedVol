// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-01-26
 *
 * Description:
 *   Optimisation helper utilities: ATM weighting, bounds diagnostics,
 *   and logging.
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

#include "Core/Types.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>

namespace uv::math::opt
{
/**
 * @brief Parameters for the ATM weighting kernel.
 *
 * @details
 * Defines the amplitude and width of a Gaussian bump applied around ATM
 * in log-moneyness coordinates. The weight used is:
 *
 *   w(x) = sqrt( 1 + (wATM - 1) * exp( -(x/k0)^2 ) )
 *
 * where x = log(K/F).
 *
 * @tparam T Floating-point type.
 */
template <std::floating_point T> struct WeightATM
{
    T wATM{1.0}; ///< ATM weight amplitude (wATM = 1 means no reweighting).
    T k0{};      ///< Width parameter controlling how quickly the bump decays.
};

/**
 * @brief Compute ATM weights on a log(K/F) grid.
 *
 * @details
 * Fills @p out with weights computed from @p logKF and @p params:
 *
 *   out[i] = sqrt( 1 + (wATM - 1) * exp( -(logKF[i]/k0)^2 ) )
 *
 * The output weight is intended to scale residuals (and their Jacobians)
 * in a weighted least squares objective.
 *
 * @tparam T Floating-point type.
 * @param logKF Log-moneyness grid values log(K/F).
 * @param params ATM kernel parameters.
 * @param out Output weights, same length as @p logKF.
 * @param doValidate Enable input validation checks.
 */
template <std::floating_point T>
void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate = true
);

/**
 * @brief Clamp an initial guess inside component-wise bounds.
 *
 * @details
 * For each component i:
 *
 *   initGuess[i] = clamp(initGuess[i], lowerBounds[i], upperBounds[i])
 *
 * Logs a warning for any parameter that is modified.
 *
 * @tparam N Number of parameters.
 * @param initGuess Initial parameter vector (modified in-place).
 * @param lowerBounds Component-wise lower bounds.
 * @param upperBounds Component-wise upper bounds.
 * @param paramNames Parameter names used for logging.
 */
template <std::size_t N>
void clamp(
    std::array<double, N>& initGuess,
    const std::array<double, N>& lowerBounds,
    const std::array<double, N>& upperBounds,
    const std::array<std::string_view, N>& paramNames
) noexcept;

/**
 * @brief Emit warnings when parameters lie numerically on their bounds.
 *
 * @details
 * Uses small absolute/relative tolerances to detect whether a value is
 * effectively equal to its lower or upper bound, and logs warnings.
 *
 * @tparam N Number of parameters.
 * @param x Optimised parameter vector.
 * @param lowerBounds Component-wise lower bounds.
 * @param upperBounds Component-wise upper bounds.
 * @param paramNames Parameter names used for logging.
 */
template <std::size_t N>
void warnBoundsHit(
    std::span<double> x,
    const std::array<double, N>& lowerBounds,
    const std::array<double, N>& upperBounds,
    const std::array<std::string_view, N>& paramNames
) noexcept;

/**
 * @brief Log optimisation results in a standardised calibration format.
 *
 * @details
 * Prints parameter values and summary diagnostics, including SSE, elapsed
 * time, iteration count, and a success flag.
 *
 * @tparam N Number of parameters.
 * @param x Optimised parameter vector.
 * @param paramNames Parameter names used for logging.
 * @param sse Final sum of squared errors.
 * @param iterCount Iteration count.
 * @param elapsedMs Elapsed time in milliseconds.
 * @param isSuccess Whether the optimiser reported success.
 */
template <std::size_t N>
void logResults(
    std::span<double> x,
    const std::array<std::string_view, N>& paramNames,
    double sse,
    unsigned iterCount,
    double elapsedMs,
    bool isSuccess
) noexcept;

namespace detail
{
/**
 * @brief Validate inputs for weightsATM.
 *
 * @details
 * Checks size consistency and basic parameter validity. Intended for
 * debug or user-facing validation prior to computing weights.
 *
 * @tparam T Floating-point type.
 * @param logKF Log-moneyness grid values.
 * @param params ATM kernel parameters.
 * @param out Output buffer for weights.
 */
template <std::floating_point T>
void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
);
} // namespace detail
} // namespace uv::math::opt

#include "Functions.inl"