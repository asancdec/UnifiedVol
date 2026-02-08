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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

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

template <opt::nlopt::Algorithm Algo>
void fillCalendarMContext(
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
