// SPDX-License-Identifier: Apache-2.0

#include "Budgets.hpp"
#include "Math/Interpolation/BSpline/Interpolator.hpp"
#include "Timing.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <vector>

namespace
{
std::vector<double> makeOpenUniformCubicKnots(const std::size_t controlPointCount)
{
    constexpr std::size_t degree{3};

    std::vector<double> knots;
    knots.reserve(controlPointCount + degree + 1);

    for (std::size_t i{0}; i <= degree; ++i)
    {
        knots.emplace_back(0.0);
    }

    const std::size_t interiorCount{controlPointCount - degree - 1};
    for (std::size_t i{1}; i <= interiorCount; ++i)
    {
        knots.emplace_back(static_cast<double>(i));
    }

    const double rightEndpoint{static_cast<double>(interiorCount + 1)};
    for (std::size_t i{0}; i <= degree; ++i)
    {
        knots.emplace_back(rightEndpoint);
    }

    return knots;
}

std::vector<double> makeControlPoints(const std::size_t controlPointCount)
{
    std::vector<double> controlPoints;
    controlPoints.reserve(controlPointCount);

    for (std::size_t i{0}; i < controlPointCount; ++i)
    {
        const double x{static_cast<double>(i)};
        controlPoints.emplace_back(
            0.20 + 0.05 * std::sin(0.17 * x) + 0.03 * std::cos(0.11 * x)
        );
    }

    return controlPoints;
}

std::vector<double>
makeGrid(const double left, const double right, const std::size_t size)
{
    std::vector<double> x(size);

    const double dx{(right - left) / static_cast<double>(size - 1)};
    for (std::size_t i{0}; i < size; ++i)
    {
        x[i] = left + static_cast<double>(i) * dx;
    }

    return x;
}

} // namespace

TEST(PerformanceBSpline, EvaluatesLargeCubicGridWithinThroughputBudget)
{
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        "bsplineLargeEvaluation"
    );

    constexpr std::size_t controlPointCount{64};
    constexpr std::size_t gridSize{200'000};

    const auto controlPoints{makeControlPoints(controlPointCount)};
    const auto knots{makeOpenUniformCubicKnots(controlPointCount)};

    const double left{knots[3]};
    const double right{knots[controlPointCount]};

    const auto x{makeGrid(left, right, gridSize)};

    const uv::math::interp::bspline::BSpline<double, 3> spline{controlPoints, knots};

    std::vector<double> out(x.size());

    spline.evalInplace(out, x);

    ASSERT_EQ(out.size(), x.size());
    EXPECT_TRUE(std::isfinite(out.front()));
    EXPECT_TRUE(std::isfinite(out[out.size() / 2]));
    EXPECT_TRUE(std::isfinite(out.back()));

    const double ms = uv::tests::performance::bestElapsedMs(
        [&]
        {
            spline.evalInplace(out, x);
        }
    );

    EXPECT_LT(ms, budget.maxMs);
}
