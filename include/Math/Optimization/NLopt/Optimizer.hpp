
// SPDX-License-Identifier: Apache-2.0
/*
 * File:        double.hpp
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

#include "Math/Optimization/NLopt/Config.hpp"
#include "Utils/Aux/StopWatch.hpp"

#include <nlopt.hpp>
#include <array>
#include <optional>
#include <cstddef>

namespace uv::math::opt::nlopt
{
    /**
     * @brief Lightweight RAII wrapper around an NLopt optimizer.
     *
     * This class provides a structured interface for configuring and running
     * NLopt optimizations with:
     *
     * - fixed parameter dimension at compile time
     * - bound handling and clamping of initial guesses
     * - inequality constraints with shared tolerance
     * - analytical objective gradients
     * - iteration counting and timing
     * - optional user-defined auxiliary scalar value
     *
     * The optimizer is configured once via @ref Config and can be reused
     * by spawning fresh instances using @ref fresh().
     *
     * @tparam N     Number of optimization parameters.
     * @tparam Algo  NLopt algorithm identifier.
     *
     * @note
     * - NLopt callbacks require non-const `void*` user data; this class ensures
     *   that all referenced contexts remain alive for the duration of the solve.
     * - This class is not thread-safe.
     */
    template <std::size_t N, ::nlopt::algorithm Algo>
    class Optimizer
    {
    private:
        /// NLopt-compatible objective / constraint function signature
        using NloptFunction = double (*)(unsigned, const double*, double*, void*);

        //--------------------------------------------------------------------------
        // Configuration and engine
        //--------------------------------------------------------------------------

        Config<N>        config_;   ///< Numerical tolerances and metadata
        ::nlopt::opt     opt_;      ///< Underlying NLopt optimizer
        utils::StopWatch timer_;    ///< Execution timer

        //--------------------------------------------------------------------------
        // Bounds and initial guess
        //--------------------------------------------------------------------------

        std::array<double, N> lowerBounds_; ///< Lower parameter bounds
        std::array<double, N> upperBounds_; ///< Upper parameter bounds
        std::array<double, N> initGuess_;   ///< Initial guess (clamped to bounds)

        //--------------------------------------------------------------------------
        // Objective state
        //--------------------------------------------------------------------------

        NloptFunction userFn_;       ///< User-provided objective function
        void* userData_;             ///< User-provided context pointer
        unsigned      iterCount_;    ///< Number of objective evaluations

        //--------------------------------------------------------------------------
        // Generic auxiliary storage
        //--------------------------------------------------------------------------

        std::optional<double> userValue_; ///< Optional scalar value shared with callbacks

        /**
         * @brief Static thunk routing NLopt callbacks to this instance.
         *
         * This function increments the iteration counter and forwards
         * evaluation to the user-provided objective function.
         */
        static double ObjectiveThunk(
            unsigned n,
            const double* x,
            double* grad,
            void* p
        ) noexcept;

    public:
        //--------------------------------------------------------------------------
        // Construction
        //--------------------------------------------------------------------------

        /// Deleted default constructor
        Optimizer() = delete;

        /**
         * @brief Construct an optimizer from a configuration object.
         *
         * Sets global NLopt parameters such as maximum evaluations and
         * relative objective tolerance.
         *
         * @param config Optimizer configuration.
         */
        explicit Optimizer(const Config<N>& config);

        /**
         * @brief Create a fresh optimizer with identical configuration.
         *
         * Returns a new optimizer instance sharing the same configuration
         * but with no constraints, objective, or state set.
         *
         * @return New optimizer instance.
         */
        Optimizer<N, Algo> fresh() const noexcept;

        //--------------------------------------------------------------------------
        // Setup
        //--------------------------------------------------------------------------

        /**
         * @brief Set initial guess and parameter bounds.
         *
         * The initial guess is clamped to the provided bounds before being
         * passed to NLopt.
         *
         * @param initGuess   Initial parameter guess.
         * @param lowerBounds Lower bounds for parameters.
         * @param upperBounds Upper bounds for parameters.
         */
        void setGuessBounds(
            std::array<double, N> initGuess,
            std::array<double, N> lowerBounds,
            std::array<double, N> upperBounds
        ) noexcept;

        /**
         * @brief Add an inequality constraint.
         *
         * Adds a constraint of the form c(x) <= 0 using the configured
         * constraint tolerance.
         *
         * @param c     Constraint function.
         * @param data  User-defined context (must remain alive during optimization).
         */
        void addInequalityConstraint(NloptFunction c, void* data) noexcept;

        /**
         * @brief Set the objective function to be minimized.
         *
         * Resets the iteration counter and registers the user-provided
         * objective function and context.
         *
         * @param f     Objective function.
         * @param data  User-defined context (must remain alive during optimization).
         */
        void setMinObjective(NloptFunction f, void* data) noexcept;

        //--------------------------------------------------------------------------
        // Execution
        //--------------------------------------------------------------------------

        /**
         * @brief Run the optimization.
         *
         * Executes the NLopt solver using the configured bounds, constraints,
         * and objective function.
         *
         * Logs results, timing, and boundary hits before returning.
         *
         * @return Optimized parameter vector.
         */
        Vector<double> optimize();

        //--------------------------------------------------------------------------
        // Auxiliary user value
        //--------------------------------------------------------------------------

        /**
         * @brief Set a user-defined scalar value.
         *
         * This value can be accessed by objective or constraint callbacks
         * via @ref userValue().
         *
         * @param v Value to store.
         */
        void setUserValue(double v) noexcept;

        //--------------------------------------------------------------------------
        // Accessors
        //--------------------------------------------------------------------------

        /// @return Inequality epsilon used in constraints.
        const double& eps() const noexcept;

        /// @return Constraint tolerance.
        double tol() const noexcept;

        /**
         * @brief Retrieve the user-defined scalar value.
         *
         * @throws ErrorCode::InvalidArgument if the value was not set.
         */
        const double& userValue() const noexcept;
    };
}  // namespace uv::math::opt::nlopt 


#include "Optimizer.inl"