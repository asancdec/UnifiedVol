// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-01-26
 *
 * Description:
 *   Optimization helper utilities.
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

#include <array>
#include <concepts>
#include <cstddef>
#include <vector>
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
 * @brief Clamp parameters to their bounds.
 *
 * Optionally validates input sizes and bounds, then clamps each value
 * in @p initGuess to the corresponding [lowerBounds, upperBounds] interval.
 *
 * @param initGuess     Initial parameter guesses (modified in place).
 * @param lowerBounds   Lower bounds per parameter.
 * @param upperBounds   Upper bounds per parameter.
 * @param doValidate    Enable input validation.
 */
void clamp(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate = true
);

/**
 * @brief Emit warnings when parameters are close to their bounds.
 *
 * Checks whether any value in @p x is numerically close to its lower or
 * upper bound and logs a warning if so.
 *
 * @param x            Parameter values to inspect.
 * @param lowerBounds  Lower bounds per parameter.
 * @param upperBounds  Upper bounds per parameter.
 * @param doValidate   Enable input validation.
 */
void warnBoundsHit(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate = true
);

/**
 * @brief Log calibration results.
 *
 * Logs parameter values (if names are provided), final error metric,
 * timing information, iteration count, and success status.
 *
 * @param x           Final parameter values.
 * @param paramNames  Optional parameter names (must match @p x if non-empty).
 * @param sse         Final sum of squared errors.
 * @param iterCount   Number of iterations performed.
 * @param elapsedMs   Elapsed wall-clock time in milliseconds.
 * @param isSuccess   Whether calibration converged successfully.
 */
void logResults(
    std::span<const double> x,
    std::span<const std::string_view> paramNames,
    double sse,
    unsigned iterCount,
    double elapsedMs,
    bool isSuccess
);

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

void validateBounds(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
);


} // namespace detail
} // namespace uv::math::opt

#include "Functions.inl"