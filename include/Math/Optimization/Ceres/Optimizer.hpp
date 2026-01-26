// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Optimizer.hpp
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

#include "Math/Optimization/Ceres/Config.hpp"
#include "Policy.hpp"

#include <array>
#include <ceres/ceres.h>

namespace uv::math::opt::ceres
{
/**
 * @brief Lightweight wrapper around a Ceres least-squares optimisation problem.
 *
 * This class provides a structured interface for configuring and running
 * Ceres-based calibrations with:
 *
 * - fixed parameter dimension at compile time
 * - bound handling and clamping of initial guesses
 * - analytic residual blocks (user-provided CostFunction)
 * - optional robust loss via @ref Policy
 * - standardised logging of results and boundary hits
 *
 * The optimiser is configured once via @ref Config.
 *
 * @tparam N       Number of optimisation parameters.
 * @tparam Policy  Compile-time policy selecting solver strategy, linear solver,
 *                 and robust-loss construction (see @ref Policy).
 *
 * @note
 * - This class is not thread-safe.
 * - Bounds are enforced via Ceres parameter bounds on a single parameter block.
 */
template <std::size_t N, typename Policy = Policy<>> class Optimizer
{
  private:
    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------

    Config config_;                      // Optimiser configuration
    std::array<double, N> lowerBounds_; // Lower parameter bounds
    std::array<double, N> upperBounds_; // Upper parameter bounds
    std::array<double, N> x_;           // Parameter block (optimised in-place)
    ::ceres::Problem problem_;          // Ceres problem instance

  public:
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    /// Deleted default constructor.
    Optimizer() = delete;

    /**
     * @brief Construct an optimiser from a configuration object.
     *
     * @param config Optimiser configuration (tolerances, parameter names,
     * logging options).
     */
    explicit Optimizer(const Config& config);

    //--------------------------------------------------------------------------
    // Calibration
    //--------------------------------------------------------------------------

    /**
     * @brief Set the initial guess and per-parameter bounds.
     *
     * Stores the input vectors, clamps the initial guess to the bounds,
     * registers a single parameter block with Ceres, and applies parameter
     * bounds on that same block.
     *
     * @param initGuess   Initial parameter guess.
     * @param lowerBounds Lower parameter bounds.
     * @param upperBounds Upper parameter bounds.
     *
     * @note The initial guess is clamped in-place before being passed to Ceres.
     */
    void setGuessBounds(
        const std::array<double, N>& initGuess,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds
    ) noexcept;

    /**
     * @brief Add an analytic residual block to the problem.
     *
     * Adds a user-provided Ceres cost function (analytic Jacobian) to the
     * internal problem, optionally wrapped by a robust loss function created
     * by the @ref Policy.
     *
     * @param cf Owning pointer to a Ceres cost function.
     *
     * @note Ownership of the cost function is transferred to Ceres.
     * @note The robust loss (if any) is constructed using @c config_.lossScale.
     */
    void addAnalyticResidual(std::unique_ptr<::ceres::CostFunction> cf) noexcept;

    /**
     * @brief Solve the optimisation problem.
     *
     * Builds @c ceres::Solver::Options from @ref Config and @ref Policy,
     * executes @c ceres::Solve, then:
     *
     * - warns if parameters are (numerically) on bounds
     * - logs the final calibration line (params, SSE, time, iterations, status)
     *
     * @return Optimised parameter vector.
     */
    std::array<double, N> optimize();
};
} // namespace uv::math::opt::ceres

#include "Optimizer.inl"