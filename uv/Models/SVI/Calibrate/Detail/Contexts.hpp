// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Models/SVI/Calibrate/Detail/SliceData.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <cstddef>
#include <span>

namespace uv::models::svi::detail
{

struct ObjectiveContexts
{
    const double* k;
    const double* wM;
    std::size_t n;
    const double atmTotalVariance;

    explicit ObjectiveContexts(
        std::span<const double> logKF,
        std::span<const double> totalVariance,
        double atmTotalVariance
    ) noexcept;
};

struct CalendarMContext
{
    std::span<const double> logKF;
    std::span<const double> prevWk;
    double atmTotalVariance;
    double eps;
};

struct SliceConstraints
{
    CalendarMContext calM;
    Vector<double> calLogKF;
    Vector<double> calPrevWk;
};

struct ConvexityMContext
{
    std::span<const double> logKF;
    double atmTotalVariance{};
};

template <opt::nlopt::Algorithm Algo> void fillCalendarMContext(
    CalendarMContext& mctx,
    Vector<double>& logKFBuf,
    Vector<double>& prevWkBuf,
    std::size_t numStrikes,
    const Params<double>& prevParams,
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    std::span<const double> logKF,
    const SliceData& sliceData
) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Contexts.inl"
