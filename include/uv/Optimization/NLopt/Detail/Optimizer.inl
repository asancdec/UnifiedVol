// SPDX-License-Identifier: Apache-2.0
/*
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

#include "Base/Errors/Errors.hpp"
#include "Base/Macros/Require.hpp"
#include "Optimization/Helpers.hpp"
#include "Optimization/NLopt/Detail/MapAlgorithm.hpp"
#include "Optimization/NLopt/Detail/NLoptStatus.hpp"
#include <nlopt.hpp>
#include <string_view>

#include <exception> // std::exception
#include <iostream>
#include <iostream> // std::cerr
#include <stdexcept>

namespace uv::opt::nlopt
{
template <std::size_t N, Algorithm Algo>
Optimizer<N, Algo>::Optimizer(const Config<N>& config)
    : config_(config),
      opt_(detail::toNlopt(Algo), N),

      lowerBounds_(),
      upperBounds_(),
      initGuess_(),
      userFn_(nullptr),
      userData_(nullptr),
      iterCount_(0U)
{
    opt_.set_ftol_rel(config_.ftolRel);
    opt_.set_maxeval(config_.maxEval);
    opt_.set_exceptions_enabled(false);
}

template <std::size_t N, Algorithm Algo>
Optimizer<N, Algo> Optimizer<N, Algo>::fresh() const noexcept
{
    return Optimizer<N, Algo>{config_};
}

template <std::size_t N, Algorithm Algo> void Optimizer<N, Algo>::setGuessBounds(
    std::array<double, N> initGuess,
    std::array<double, N> lowerBounds,
    std::array<double, N> upperBounds
) noexcept
{

    clampBounds(initGuess, lowerBounds, upperBounds);

    initGuess_ = initGuess;
    lowerBounds_ = lowerBounds;
    upperBounds_ = upperBounds;

    opt_.set_lower_bounds(Vector<double>(lowerBounds_.begin(), lowerBounds_.end()));
    opt_.set_upper_bounds(Vector<double>(upperBounds_.begin(), upperBounds_.end()));
}

template <std::size_t N, Algorithm Algo> double Optimizer<N, Algo>::objectiveThunk(
    unsigned n,
    const double* x,
    double* grad,
    void* p
) noexcept
{
    auto* self = static_cast<Optimizer<N, Algo>*>(p);
    ++self->iterCount_;
    return self->userFn_ ? self->userFn_(n, x, grad, self->userData_) : 0.0;
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::addInequalityConstraint(NloptFunction c, void* data) noexcept
{
    opt_.add_inequality_constraint(c, data, config_.tol);
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::addInequalityMConstraint(
    std::size_t m,
    NloptMFunction c,
    void* data
) noexcept
{
    Vector<double> tol(m, config_.tol);
    opt_.add_inequality_mconstraint(c, data, tol);
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::setMinObjective(NloptFunction f, void* data) noexcept
{
    iterCount_ = 0U;
    userFn_ = f;
    userData_ = data;

    if (config_.verbose)
    {
        opt_.set_min_objective(&Optimizer<N, Algo>::objectiveThunk, this);
    }
    else
    {
        opt_.set_min_objective(userFn_, userData_);
    }
}

template <std::size_t N, Algorithm Algo> Vector<double> Optimizer<N, Algo>::optimize()
{

    Vector<double> x(initGuess_.cbegin(), initGuess_.cend());
    double sse{0.0};

    timer_.StartStopWatch();

    ::nlopt::result successCode = opt_.optimize(x, sse);

    timer_.StopStopWatch();

    if (config_.verbose)
    {
        warnBoundsHit(x, lowerBounds_, upperBounds_);

        logResults(
            x,
            config_.paramNames,
            sse,
            iterCount_,
            timer_.GetTime<std::milli>(),
            (successCode > ::nlopt::FAILURE),
            detail::toString(successCode)
        );
    }

    return x;
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::setUserValue(double v) noexcept
{
    userValue_ = v;
}

template <std::size_t N, Algorithm Algo>
const double& Optimizer<N, Algo>::eps() const noexcept
{
    return config_.eps;
}

template <std::size_t N, Algorithm Algo> double Optimizer<N, Algo>::tol() const noexcept
{
    return config_.tol;
}

template <std::size_t N, Algorithm Algo>
const double& Optimizer<N, Algo>::userValue() const noexcept
{
    UV_REQUIRE_SET(userValue_);

    return *userValue_;
}
} // namespace uv::opt::nlopt