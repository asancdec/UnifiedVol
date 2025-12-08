// SPDX-License-Identifier: Apache-2.0
/*
 * File:        double.inl
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


#include "Math/Optimization/Functions.hpp"

namespace uv::math::opt::nlopt
{
    template <std::size_t N, ::nlopt::algorithm Algo>
    Optimizer<N, Algo>::Optimizer(const Config<N>& config) :
        config_(config), 
        opt_(Algo, N), 
        timer_(),
        lowerBounds_(),
        upperBounds_(), 
        initGuess_(),
        userFn_(nullptr), 
        userData_(nullptr), 
        iterCount_(0U) 
    {
        opt_.set_ftol_rel(config_.ftolRel); 
        opt_.set_maxeval(config_.maxEval);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    Optimizer<N, Algo> Optimizer<N, Algo>::fresh() const noexcept
    {
        return Optimizer<N, Algo>{ config_ };
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Optimizer<N, Algo>::setGuessBounds(std::array<double, N> initGuess,
        std::array<double, N> lowerBounds,
        std::array<double, N> upperBounds) noexcept
    {
        // Clamp initial guess within upper and lower bounds
        clamp<N>(initGuess, lowerBounds, upperBounds, config_.paramNames);

        // Store arrays
        initGuess_ = initGuess;
        lowerBounds_ = lowerBounds;
        upperBounds_ = upperBounds;

        // Configure NLopt with bounds and tolerances
        opt_.set_lower_bounds(Vector<double>(lowerBounds_.begin(), lowerBounds_.end()));
        opt_.set_upper_bounds(Vector<double>(upperBounds_.begin(), upperBounds_.end()));
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    double Optimizer<N, Algo>::ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept
    {
        auto* self = static_cast<Optimizer<N, Algo>*>(p);
        ++self->iterCount_;
        return self->userFn_ ? self->userFn_(n, x, grad, self->userData_) : 0.0;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Optimizer<N, Algo>::addInequalityConstraint(
        NloptFunction c,
        void* data) noexcept
    {
        opt_.add_inequality_constraint(c, data, config_.tol);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Optimizer<N, Algo>::setMinObjective(NloptFunction f, void* data) noexcept
    {
        iterCount_ = 0U;
        userFn_ = f;
        userData_ = data;

        // Route NLopt callback into this instance
        opt_.set_min_objective(&Optimizer<N, Algo>::ObjectiveThunk, this);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    Vector<double> Optimizer<N, Algo>::optimize() noexcept
    {
        // Copy initial guess to working vector
        Vector<double> x(initGuess_.cbegin(), initGuess_.cend());
        double sse{ 0.0 };

        // Start timer
        timer_.StartStopWatch();

        // Run optimizer
        ::nlopt::result successCode = opt_.optimize(x, sse);

        // End timer
        timer_.StopStopWatch();

        // Warn if upper or lower bounds are touched
        warnBoundsHit
        (
            x,
            lowerBounds_,
            upperBounds_,
            config_.paramNames
        );

        // Log calibration results 
        logResults(
            x,                                       // Parameters
            config_.paramNames,                      // Parameter names
            sse,                                     // SSE
            iterCount_,                              // Iterations
            timer_.GetTime<std::milli>(),            // Elapsed [ms]
            (successCode > ::nlopt::FAILURE)         // Success flag
        );

        return x;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    const double& Optimizer<N, Algo>::eps() const noexcept
    {
        return config_.eps;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    double Optimizer<N, Algo>::tol() const noexcept
    {
        return config_.tol;
    }
}
