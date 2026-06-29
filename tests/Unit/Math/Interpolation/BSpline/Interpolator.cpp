// SPDX-License-Identifier: Apache-2.0

#include "Math/Interpolation/BSpline/Interpolator.hpp"
#include "Base/Errors/Errors.hpp"
#include "Support/Tolerances.hpp"

#include <algorithm>
#include <cstddef>
#include <gtest/gtest.h>
#include <vector>

namespace bspline = uv::math::interp::bspline;

TEST(MathBSplineInterpolator, LinearSplineInterpolatesControlPolygon)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    const bspline::BSpline<double, 1> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.5, 1.0, 1.5, 2.0};

    const auto y{spline.eval(x)};

    ASSERT_EQ(y.size(), x.size());
    EXPECT_NEAR(y[0], 0.0, 1e-15);
    EXPECT_NEAR(y[1], 1.0, 1e-15);
    EXPECT_NEAR(y[2], 2.0, 1e-15);
    EXPECT_NEAR(y[3], 3.0, 1e-15);
    EXPECT_NEAR(y[4], 4.0, 1e-15);
}

TEST(MathBSplineInterpolator, UnitControlSplinesFormPartitionOfUnityOnClampedDomain)
{
    constexpr std::size_t degree{3};

    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    const std::vector<double> x{0.0, 0.2, 0.7, 1.0, 1.5, 2.25, 2.9, 3.0};

    constexpr std::size_t basisCount{6};

    std::vector<double> sum(x.size(), 0.0);

    for (std::size_t basisIdx{0}; basisIdx < basisCount; ++basisIdx)
    {
        std::vector<double> controlPoints(basisCount, 0.0);
        controlPoints[basisIdx] = 1.0;

        std::vector<double> basis(x.size());

        const bspline::BSpline<double, degree> basisSpline{controlPoints, knots};

        basisSpline.evalInplace(basis, x);

        for (std::size_t i{0}; i < x.size(); ++i)
        {
            EXPECT_GE(basis[i], -1e-15) << "basisIdx=" << basisIdx << ", x=" << x[i];
            EXPECT_LE(basis[i], 1.0 + 1e-15) << "basisIdx=" << basisIdx << ", x=" << x[i];

            sum[i] += basis[i];
        }
    }

    for (std::size_t i{0}; i < x.size(); ++i)
    {
        EXPECT_NEAR(sum[i], 1.0, 1e-14) << "x=" << x[i];
    }
}

TEST(MathBSplineInterpolator, ConstantControlPointsEvaluateToConstant)
{
    const std::vector<double> controlPoints(6, 7.5);
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    const bspline::BSpline<double, 3> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.25, 0.75, 1.5, 2.5, 3.0};

    const auto y{spline.eval(x)};

    for (double value : y)
    {
        EXPECT_NEAR(value, 7.5, 1e-14);
    }
}

TEST(MathBSplineInterpolator, ClampedSplineStaysWithinControlPointBounds)
{
    const std::vector<double> controlPoints{-1.0, 0.5, 3.0, 2.0, 4.5, 1.5};
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    const bspline::BSpline<double, 3> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.1, 0.4, 0.9, 1.4, 2.0, 2.6, 3.0};

    const auto [minIt, maxIt] =
        std::minmax_element(controlPoints.begin(), controlPoints.end());

    const auto y{spline.eval(x)};

    for (std::size_t i{0}; i < y.size(); ++i)
    {
        EXPECT_GE(y[i] + uv::tests::tolerance::InterpolationInvariant, *minIt)
            << "x=" << x[i];

        EXPECT_LE(y[i], *maxIt + uv::tests::tolerance::InterpolationInvariant)
            << "x=" << x[i];
    }
}

TEST(MathBSplineInterpolator, EvalInplaceMatchesEval)
{
    const std::vector<double> controlPoints{0.0, 1.0, 0.0, 2.0};
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0};
    const bspline::BSpline<double, 2> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.5, 1.0, 1.5, 2.0};

    std::vector<double> out(x.size());

    spline.evalInplace(out, x);
    const auto y{spline.eval(x)};

    EXPECT_EQ(out, y);
}

TEST(MathBSplineInterpolator, SetControlPointsUpdatesOwnedState)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    bspline::BSpline<double, 1> spline{controlPoints, knots};
    std::vector<double> nextControlPoints{1.0, 3.0, 5.0};

    spline.setControlPoints(nextControlPoints);
    nextControlPoints.assign(nextControlPoints.size(), 100.0);

    const auto y{spline.eval(std::vector<double>{0.0, 1.0, 2.0})};

    ASSERT_EQ(y.size(), 3U);
    EXPECT_DOUBLE_EQ(y[0], 1.0);
    EXPECT_DOUBLE_EQ(y[1], 3.0);
    EXPECT_DOUBLE_EQ(y[2], 5.0);
}

TEST(MathBSplineInterpolator, SetKnotsUpdatesOwnedState)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    bspline::BSpline<double, 1> spline{controlPoints, knots};
    std::vector<double> nextKnots{0.0, 0.0, 2.0, 4.0, 4.0};

    spline.setKnots(nextKnots);
    nextKnots.assign(nextKnots.size(), -10.0);

    const auto y{spline.eval(std::vector<double>{1.0, 2.0, 3.0})};

    ASSERT_EQ(y.size(), 3U);
    EXPECT_DOUBLE_EQ(y[0], 1.0);
    EXPECT_DOUBLE_EQ(y[1], 2.0);
    EXPECT_DOUBLE_EQ(y[2], 3.0);
}

TEST(MathBSplineInterpolator, RejectedSettersPreservePreviousState)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    bspline::BSpline<double, 1> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 1.0, 2.0};
    const auto before{spline.eval(x)};

    EXPECT_THROW(
        spline.setControlPoints(std::vector<double>{1.0, 2.0}),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        spline.setKnots(std::vector<double>{0.0, 0.0, 2.0, 1.0, 2.0}),
        uv::errors::UnifiedVolError
    );

    EXPECT_EQ(spline.eval(x), before);
}

TEST(MathBSplineInterpolator, EvalInplaceMatchesIndependentScalarEvaluations)
{
    const std::vector<double> controlPoints{-1.0, 0.5, 3.0, 2.0, 4.5, 1.5};
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    const bspline::BSpline<double, 3> spline{controlPoints, knots};

    const std::vector<double>
        x{0.0, 0.1, 1.0, 1.0, 0.4, -0.2, 0.9, 1.4, 3.1, 2.0, 2.6, 3.0};
    std::vector<double> out(x.size());

    spline.evalInplace(out, x);

    for (std::size_t i{0}; i < x.size(); ++i)
    {
        const auto value{spline.eval(std::vector<double>{x[i]})};

        ASSERT_EQ(value.size(), 1U);
        EXPECT_NEAR(value.front(), out[i], 1e-15) << "x=" << x[i];
    }
}

TEST(MathBSplineInterpolator, EvalInplaceMatchesEvalForSortedInputs)
{
    const std::vector<double> controlPoints{-1.0, 0.5, 3.0, 2.0, 4.5, 1.5};
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    const bspline::BSpline<double, 3> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.1, 0.4, 0.9, 1.4, 2.0, 2.6, 3.0};
    std::vector<double> inplaceEval(x.size());

    spline.evalInplace(inplaceEval, x);
    const auto vectorEval{spline.eval(x)};

    ASSERT_EQ(inplaceEval.size(), vectorEval.size());

    for (std::size_t i{0}; i < x.size(); ++i)
    {
        EXPECT_NEAR(inplaceEval[i], vectorEval[i], 1e-15) << "x=" << x[i];
    }
}

TEST(MathBSplineInterpolator, EvalReturnsZeroOutsideDomain)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    const bspline::BSpline<double, 1> spline{controlPoints, knots};

    const auto y{spline.eval(std::vector<double>{-0.1, 2.1})};

    ASSERT_EQ(y.size(), 2U);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(y[1], 0.0);
}

TEST(MathBSplineInterpolator, EvalAtRightEndpointReturnsLastControlPoint)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    const bspline::BSpline<double, 1> spline{controlPoints, knots};

    const auto y{spline.eval(std::vector<double>{2.0})};

    ASSERT_EQ(y.size(), 1U);
    EXPECT_DOUBLE_EQ(y.front(), 4.0);
}

TEST(MathBSplineInterpolator, DegreeZeroSplineIsPiecewiseConstant)
{
    const std::vector<double> controlPoints{1.0, 3.0, 5.0};
    const std::vector<double> knots{0.0, 1.0, 2.0, 3.0};
    const bspline::BSpline<double, 0> spline{controlPoints, knots};
    const std::vector<double> x{0.25, 1.25, 2.25, 3.0};

    const auto y{spline.eval(x)};

    ASSERT_EQ(y.size(), x.size());
    EXPECT_DOUBLE_EQ(y[0], 1.0);
    EXPECT_DOUBLE_EQ(y[1], 3.0);
    EXPECT_DOUBLE_EQ(y[2], 5.0);
    EXPECT_DOUBLE_EQ(y[3], 5.0);
}

TEST(MathBSplineInterpolator, RejectsInconsistentKnotCount)
{
    const std::vector<double> controlPoints{1.0, 2.0, 3.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};

    EXPECT_THROW(
        (bspline::BSpline<double, 2>{controlPoints, knots}),
        uv::errors::UnifiedVolError
    );
}

TEST(MathBSplineInterpolator, RejectsNonDecreasingViolationsAndDegenerateDomain)
{
    const std::vector<double> controlPoints{1.0, 2.0, 3.0};

    EXPECT_THROW(
        (bspline::BSpline<double, 1>{
            controlPoints,
            std::vector<double>{0.0, 0.0, 2.0, 1.0, 2.0}
        }),
        uv::errors::UnifiedVolError
    );

    EXPECT_THROW(
        (bspline::BSpline<double, 1>{
            controlPoints,
            std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0}
        }),
        uv::errors::UnifiedVolError
    );
}

TEST(MathBSplineInterpolator, RejectsEvalInplaceSizeMismatch)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    const bspline::BSpline<double, 1> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> out(2);

    EXPECT_THROW(spline.evalInplace(out, x), uv::errors::UnifiedVolError);
}
