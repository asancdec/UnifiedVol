// SPDX-License-Identifier: Apache-2.0

#include "Math/Interpolation/BSpline/Interpolator.hpp"
#include "Base/Errors/Errors.hpp"
#include "Math/Interpolation/BSpline/Detail/Evaluate.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(MathBSplineInterpolator, LinearSplineInterpolatesControlPolygon)
{
    const std::vector<double> controlPoints{0.0, 2.0, 4.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 2.0, 2.0};
    const uv::math::interp::bspline::BSpline<double, 1> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.5, 1.0, 1.5, 2.0};

    const auto y = spline.eval(x);

    ASSERT_EQ(y.size(), x.size());
    EXPECT_NEAR(y[0], 0.0, 1e-15);
    EXPECT_NEAR(y[1], 1.0, 1e-15);
    EXPECT_NEAR(y[2], 2.0, 1e-15);
    EXPECT_NEAR(y[3], 3.0, 1e-15);
    EXPECT_NEAR(y[4], 4.0, 1e-15);
}

TEST(MathBSplineInterpolator, BasisFunctionsFormPartitionOfUnityOnClampedDomain)
{
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    const std::vector<double> x{0.0, 0.2, 0.7, 1.0, 1.5, 2.25, 2.9, 3.0};

    for (std::size_t basisIdx = 0; basisIdx < 6; ++basisIdx)
    {
        std::vector<double> basis(x.size());
        uv::math::interp::bspline::detail::coxDeBoor<double, 3>(
            basis,
            x,
            knots,
            basisIdx
        );

        for (double value : basis)
        {
            EXPECT_GE(value, -1e-15);
            EXPECT_LE(value, 1.0 + 1e-15);
        }
    }

    for (std::size_t sampleIdx = 0; sampleIdx < x.size() - 1; ++sampleIdx)
    {
        double sum = 0.0;
        for (std::size_t basisIdx = 0; basisIdx < 6; ++basisIdx)
        {
            std::vector<double> basis(x.size());
            uv::math::interp::bspline::detail::coxDeBoor<double, 3>(
                basis,
                x,
                knots,
                basisIdx
            );
            sum += basis[sampleIdx];
        }

        EXPECT_NEAR(sum, 1.0, 1e-14);
    }
}

TEST(MathBSplineInterpolator, ConstantControlPointsEvaluateToConstant)
{
    const std::vector<double> controlPoints(6, 7.5);
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    const uv::math::interp::bspline::BSpline<double, 3> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.25, 0.75, 1.5, 2.5, 3.0};

    const auto y = spline.eval(x);

    for (double value : y)
    {
        EXPECT_NEAR(value, 7.5, 1e-14);
    }
}

TEST(MathBSplineInterpolator, EvalInplaceMatchesEval)
{
    const std::vector<double> controlPoints{0.0, 1.0, 0.0, 2.0};
    const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0};
    const uv::math::interp::bspline::BSpline<double, 2> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 0.5, 1.0, 1.5, 2.0};
    std::vector<double> out(x.size());

    spline.evalInplace(out, x);
    const auto y = spline.eval(x);

    EXPECT_EQ(out, y);
}

TEST(MathBSplineInterpolator, RejectsInconsistentKnotCount)
{
    const std::vector<double> controlPoints{1.0, 2.0, 3.0};
    const std::vector<double> knots{0.0, 0.0, 1.0, 1.0};

    EXPECT_THROW(
        (uv::math::interp::bspline::BSpline<double, 2>{controlPoints, knots}),
        uv::errors::UnifiedVolError
    );
}

TEST(MathBSplineInterpolator, RejectsNonDecreasingViolationsAndDegenerateDomain)
{
    const std::vector<double> controlPoints{1.0, 2.0, 3.0};

    EXPECT_THROW(
        (uv::math::interp::bspline::BSpline<double, 1>{
            controlPoints,
            std::vector<double>{0.0, 0.0, 2.0, 1.0, 2.0}
        }),
        uv::errors::UnifiedVolError
    );

    EXPECT_THROW(
        (uv::math::interp::bspline::BSpline<double, 1>{
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
    const uv::math::interp::bspline::BSpline<double, 1> spline{controlPoints, knots};
    const std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> out(2);

    EXPECT_THROW(spline.evalInplace(out, x), uv::errors::UnifiedVolError);
}

TEST(MathBSplineInterpolator, DegreeZeroSplineIsPiecewiseConstant)
{
    const std::vector<double> controlPoints{1.0, 3.0, 5.0};
    const std::vector<double> knots{0.0, 1.0, 2.0, 3.0};
    const uv::math::interp::bspline::BSpline<double, 0> spline{controlPoints, knots};
    const std::vector<double> x{0.25, 1.25, 2.25, 3.0};

    const auto y = spline.eval(x);

    ASSERT_EQ(y.size(), x.size());
    EXPECT_DOUBLE_EQ(y[0], 1.0);
    EXPECT_DOUBLE_EQ(y[1], 3.0);
    EXPECT_DOUBLE_EQ(y[2], 5.0);
    EXPECT_DOUBLE_EQ(y[3], 5.0);
}
