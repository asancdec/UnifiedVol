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

#include "Models/SVI/Math.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace uv::models::svi::detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo>
void addCalendarConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    SliceConstraints& c,
    const Params<T>* prevParams,
    std::span<const double> logKF,
    const SliceData& sliceData
) noexcept
{
    if (!prevParams)
        return;

    const std::size_t numStrikes{logKF.size()};
    const Params<double> prevParamsD{prevParams->template as<double>()};

    fillCalendarMContext<Algo>(
        c.calM,
        c.calLogKF,
        c.calPrevWk,
        numStrikes,
        prevParamsD,
        optimizer,
        logKF,
        sliceData
    );

    optimizer.addInequalityMConstraint(c.calLogKF.size(), &calendarMConstraint, &c.calM);
}

template <opt::nlopt::Algorithm Algo>
void addConvexityConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    ConvexityMContext& convexityCtx,
    std::span<const double> logKF,
    double atmTotalVariance
) noexcept
{
    convexityCtx.logKF = logKF;
    convexityCtx.atmTotalVariance = atmTotalVariance;

    optimizer
        .addInequalityMConstraint(logKF.size(), &convexityMConstraint, &convexityCtx);
}

template <opt::nlopt::Algorithm Algo>
void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{
    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
        {
            const auto& opt{*static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data)};

            const double eps{opt.eps()};
            const double atmTotalVariance{opt.userValue()};
            const double b{x[0]};
            const double rho{x[1]};
            const double m{x[2]};
            const double sigma{x[3]};

            const double s2{sigma * sigma};
            const double R{std::sqrt(m * m + s2)};
            const double S{std::sqrt(1.0 - rho * rho)};

            const double base{-rho * m + R};
            const double a{atmTotalVariance - b * base};

            const double wMin{b * sigma * S + a};

            if (grad)
            {
                grad[0] = -(rho * m - R + sigma * S);
                grad[1] = -(b * (m - sigma * rho / S));
                grad[2] = -(b * (rho - m / R));
                grad[3] = -(b * (S - sigma / R));
            }

            return eps - wMin;
        },
        &optimizer
    );
}
template <opt::nlopt::Algorithm Algo>
void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
        {
            const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

            const double eps{opt.eps()};
            const double b{x[0]};
            const double onePlusRho{1.0 + x[1]};

            if (grad)
            {
                grad[0] = -onePlusRho;
                grad[1] = -b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return eps - b * onePlusRho;
        },
        &optimizer
    );

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept -> double
        {
            const auto& opt = *static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data);

            const double eps{opt.eps()};
            const double b{x[0]};
            const double oneMinusRho{1.0 - x[1]};

            if (grad)
            {
                grad[0] = -oneMinusRho;
                grad[1] = b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return eps - b * oneMinusRho;
        },
        &optimizer
    );
}

template <opt::nlopt::Algorithm Algo>
void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept
{

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) noexcept -> double
        {
            const double b{x[0]};
            const double rho{x[1]};

            if (grad)
            {
                grad[0] = (1.0 + rho);
                grad[1] = b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return b * (1.0 + rho) - 2.0;
        },
        nullptr
    );

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) noexcept -> double
        {
            const double b{x[0]};
            const double rho{x[1]};

            if (grad)
            {
                grad[0] = (1.0 - rho);
                grad[1] = -b;
                grad[2] = 0.0;
                grad[3] = 0.0;
            }

            return b * (1.0 - rho) - 2.0;
        },
        nullptr
    );
}
} // namespace uv::models::svi::detail