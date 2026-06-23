// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Math.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace uv::models::svi::detail
{

// NLopt mandates C callbacks with an opaque void* context; typed lambdas cannot
// replace these signatures without an equivalent ABI trampoline.

template <opt::nlopt::Algorithm Algo> void fillCalendarMContext(
    CalendarMContext& context,
    Vector<double>& logKBuffer,
    Vector<double>& previousVarianceBuffer,
    const Params<double>& previousParams,
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    std::span<const double> logKF,
    const SliceData& sliceData
) noexcept
{
    const std::size_t numStrikes{logKF.size()};
    const std::size_t constraintCount{numStrikes + 2};

    logKBuffer.resize(constraintCount);
    previousVarianceBuffer.resize(constraintCount);

    constexpr double wingOffset{0.15};

    for (std::size_t i = 0; i < numStrikes; ++i)
    {
        const double k{logKF[i]};
        logKBuffer[i] = k;
        previousVarianceBuffer[i] = totalVariance(
            previousParams.a,
            previousParams.b,
            previousParams.rho,
            previousParams.m,
            previousParams.sigma,
            k
        );
    }

    const auto appendWing = [&](std::size_t index, double k)
    {
        logKBuffer[index] = k;
        previousVarianceBuffer[index] = totalVariance(
            previousParams.a,
            previousParams.b,
            previousParams.rho,
            previousParams.m,
            previousParams.sigma,
            k
        );
    };

    appendWing(numStrikes, sliceData.logKFMin - wingOffset);
    appendWing(numStrikes + 1, sliceData.logKFMax + wingOffset);

    context.logKF = logKBuffer;
    context.prevWk = previousVarianceBuffer;
    context.eps = optimizer.eps();
    context.atmTotalVariance = sliceData.atmTotalVariance;
}

template <opt::nlopt::Algorithm Algo> [[gnu::hot]] double
wMinConstraint(unsigned, const double* x, double* grad, void* data) noexcept // NOSONAR
{
    const auto& optimizer{*static_cast<const opt::nlopt::Optimizer<4, Algo>*>(data)};

    const double eps{optimizer.eps()};
    const double atmTotalVariance{optimizer.userValue()};
    const double b{x[0]};
    const double rho{x[1]};
    const double m{x[2]};
    const double sigma{x[3]};

    const double sigmaSquared{sigma * sigma};
    const double radius{std::sqrt(m * m + sigmaSquared)};
    const double rhoScale{std::sqrt(1.0 - rho * rho)};

    const double base{-rho * m + radius};
    const double a{atmTotalVariance - b * base};
    const double minTotalVariance{b * sigma * rhoScale + a};

    if (grad)
    {
        grad[0] = -(rho * m - radius + sigma * rhoScale);
        grad[1] = -(b * (m - sigma * rho / rhoScale));
        grad[2] = -(b * (rho - m / radius));
        grad[3] = -(b * (rhoScale - sigma / radius));
    }

    return eps - minTotalVariance;
}

template <std::floating_point T, opt::nlopt::Algorithm Algo> void addCalendarConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    SliceConstraints& c,
    const Params<T>* prevParams,
    std::span<const double> logKF,
    const SliceData& sliceData
)
{
    if (!prevParams)
        return;

    const Params<double> prevParamsD{prevParams->template as<double>()};

    fillCalendarMContext<Algo>(
        c.calM,
        c.calLogKF,
        c.calPrevWk,
        prevParamsD,
        optimizer,
        logKF,
        sliceData
    );

    optimizer.addInequalityMConstraint(c.calLogKF.size(), &calendarMConstraint, &c.calM);
}

template <opt::nlopt::Algorithm Algo> void addConvexityConstraints(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    ConvexityMContext& convexityCtx,
    std::span<const double> logKF,
    double atmTotalVariance
)
{
    convexityCtx.logKF = logKF;
    convexityCtx.atmTotalVariance = atmTotalVariance;

    optimizer
        .addInequalityMConstraint(logKF.size(), &convexityMConstraint, &convexityCtx);
}

template <opt::nlopt::Algorithm Algo>
void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
{
    optimizer.addInequalityConstraint(&wMinConstraint<Algo>, &optimizer);
}
template <opt::nlopt::Algorithm Algo>
void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
{

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void* data) noexcept // NOSONAR
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
        +[](unsigned, const double* x, double* grad, void* data) noexcept // NOSONAR
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
void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer)
{

    optimizer.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) noexcept // NOSONAR
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
        +[](unsigned, const double* x, double* grad, void*) noexcept // NOSONAR
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
