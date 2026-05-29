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

#include "Base/Execution/ThreadPolicy.hpp"
#include "Base/Macros/Require.hpp"
#include "IO/ConsoleRedirect.hpp"
#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Helpers.hpp"

namespace uv::opt::ceres
{

template <typename PolicyT> Optimizer<PolicyT>::Optimizer(const Config& config)
    : config_(config),
      loss_(PolicyT::makeLoss(config_.lossScale))
{
    setOptions_();
}

template <typename PolicyT> void Optimizer<PolicyT>::setOptions_()
{
    options_.trust_region_strategy_type = PolicyT::trustRegionStrategy;
    options_.linear_solver_type = PolicyT::linearSolver;
    options_.max_num_iterations = config_.maxEval;
    options_.function_tolerance = config_.functionTol;
    options_.parameter_tolerance = config_.paramTol;
    options_.gradient_tolerance = config_.gradientTol;
    options_.num_threads = execution::requestThreads(config_.numThreads);
}

template <typename PolicyT> void Optimizer<PolicyT>::setBounds_(
    std::size_t n,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    const bool hasLB{!lowerBounds.empty()};
    const bool hasUB{!upperBounds.empty()};

    if (hasLB && hasUB)
    {
        validateBoundsSpec(n, lowerBounds, upperBounds);
        lowerBounds_.emplace(lowerBounds.begin(), lowerBounds.end());
        upperBounds_.emplace(upperBounds.begin(), upperBounds.end());
        return;
    }

    if (hasLB)
    {
        validateLowerBoundsSpec(n, lowerBounds);
        lowerBounds_.emplace(lowerBounds.begin(), lowerBounds.end());
    }
    else
    {
        lowerBounds_.reset();
    }

    if (hasUB)
    {
        validateUpperBoundsSpec(n, upperBounds);
        upperBounds_.emplace(upperBounds.begin(), upperBounds.end());
    }
    else
    {
        upperBounds_.reset();
    }
}

template <typename PolicyT> void Optimizer<PolicyT>::requireInitialized_() const
{
    UV_REQUIRE_VALID_STATE(
        isInitialized_,
        "Optimizer not initialized. Call initialize() first."
    );
}

template <typename PolicyT> void Optimizer<PolicyT>::requireRunStarted_() const
{
    UV_REQUIRE_VALID_STATE(isRunStarted_, "Run not started. Call beginRun() first.");
}

template <typename PolicyT> void Optimizer<PolicyT>::clampStoredBounds_()
{
    if (lowerBounds_ && upperBounds_)
    {
        clampBounds(x_, *lowerBounds_, *upperBounds_);
    }
    else if (lowerBounds_)
        clampLowerBounds(x_, *lowerBounds_);

    else if (upperBounds_)
        clampUpperBounds(x_, *upperBounds_);
}

template <typename PolicyT> void Optimizer<PolicyT>::initialize(
    std::span<const double> initGuess,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    x_.assign(initGuess.begin(), initGuess.end());
    setBounds_(x_.size(), lowerBounds, upperBounds);
    clampStoredBounds_();

    isInitialized_ = true;
    isRunStarted_ = false;
}

template <typename PolicyT> void Optimizer<PolicyT>::beginRun()
{
    requireInitialized_();

    clampStoredBounds_();

    problem_ = ::ceres::Problem{};

    const int n{static_cast<int>(x_.size())};

    problem_.AddParameterBlock(x_.data(), n);

    if (lowerBounds_)
        for (int i = 0; i < n; ++i)
            problem_.SetParameterLowerBound(x_.data(), i, (*lowerBounds_)[i]);

    if (upperBounds_)
        for (int i = 0; i < n; ++i)
            problem_.SetParameterUpperBound(x_.data(), i, (*upperBounds_)[i]);

    isRunStarted_ = true;
}

template <typename PolicyT>
void Optimizer<PolicyT>::addResidualBlock(std::unique_ptr<::ceres::CostFunction> cf)
{
    requireInitialized_();
    requireRunStarted_();

    UV_REQUIRE_NON_NULL(cf);

    problem_.AddResidualBlock(cf.release(), loss_.get(), x_.data());
}

template <typename PolicyT> void Optimizer<PolicyT>::solveInPlace()
{
    requireInitialized_();
    requireRunStarted_();

    const Verbosity v{config_.verbosity};

    ::ceres::Solver::Summary summary;

    if (v == Verbosity::FullReport)
    {
        io::ConsoleRedirect capture;
        options_.minimizer_progress_to_stdout = true;
        ::ceres::Solve(options_, &problem_, &summary);
        UV_INFO(summary.FullReport());
    }
    else
    {
        ::ceres::Solve(options_, &problem_, &summary);
    }

    if (v == Verbosity::None)
        return;

    warnBoundsHit(x_, lowerBounds_, upperBounds_);

    logResults(
        x_,
        config_.paramNames,
        summary.final_cost * 2.0,
        summary.iterations.size(),
        summary.total_time_in_seconds * 1000.0,
        (summary.termination_type == ::ceres::CONVERGENCE ||
         summary.termination_type == ::ceres::USER_SUCCESS)
    );
}

template <typename Policy> std::span<const double> Optimizer<Policy>::solve()
{
    solveInPlace();
    return params();
}

template <typename Policy> std::span<const double> Optimizer<Policy>::params() const
{
    requireInitialized_();
    requireRunStarted_();
    return std::span<const double>{x_};
}

} // namespace uv::opt::ceres
