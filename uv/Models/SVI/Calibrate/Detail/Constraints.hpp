// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <concepts>

namespace uv::models::svi::detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo> void addCalendarConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    SliceConstraints& c,
    const Params<T>* prevParams,
    std::span<const double> logKF,
    const SliceData& sliceData
);

template <opt::nlopt::Algorithm Algo> void addConvexityConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    ConvexityMContext& convexityCtxs,
    std::span<const double> logKF,
    double atmTotalVariance
);

template <opt::nlopt::Algorithm Algo>
void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

template <opt::nlopt::Algorithm Algo>
void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

template <opt::nlopt::Algorithm Algo>
void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer);

[[gnu::hot]] void calendarMConstraint(
    unsigned m,
    double* result,
    unsigned /*n*/,
    const double* x,
    double* grad,
    void* data
) noexcept;

[[gnu::hot]] void convexityMConstraint(
    unsigned m,
    double* result,
    unsigned /*n*/,
    const double* x,
    double* grad,
    void* data
) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Constraints.inl"
