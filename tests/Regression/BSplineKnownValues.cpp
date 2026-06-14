// SPDX-License-Identifier: Apache-2.0

#include "Golden.hpp"
#include "Math/Interpolation/BSpline/Interpolator.hpp"

#include <cstddef>
#include <gtest/gtest.h>

TEST(RegressionMathBSplineInterpolator, CubicSplineKnownSampleValues)
{
    const auto golden{uv::tests::golden::readBSplineKnownValues(
        "tests/Golden/bspline_known_values.json"
    )};

    const uv::math::interp::bspline::BSpline<double, 3> spline{
        golden.controlPoints,
        golden.knots
    };

    const auto y{spline.eval(golden.x)};

    ASSERT_EQ(y.size(), golden.y.size());

    for (std::size_t i{0}; i < y.size(); ++i)
    {
        EXPECT_NEAR(y[i], golden.y[i], golden.tolerance)
            << "i=" << i << " x=" << golden.x[i] << " actual=" << y[i]
            << " expected=" << golden.y[i];
    }
}