// SPDX-License-Identifier: Apache-2.0

#include "Math/Interpolation/BSpline/Interpolator.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(RegressionMathBSplineInterpolator, CubicSplineKnownSampleValues)
{
    const std::vector<double> controlPoints{0.0, 1.0, 0.0, 2.0, 1.0, 3.0};
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    const std::vector<double> x{0.0, 0.5, 1.5, 2.5, 3.0};
    const uv::math::interp::bspline::BSpline<double, 3> spline{controlPoints, knots};

    const auto y = spline.eval(x);

    ASSERT_EQ(y.size(), x.size());
    EXPECT_NEAR(y[0], 0.0, 1e-15);
    EXPECT_NEAR(y[1], 0.6354166666666666, 1e-15);
    EXPECT_NEAR(y[2], 1.0, 1e-15);
    EXPECT_NEAR(y[3], 1.4895833333333333, 1e-15);
    EXPECT_NEAR(y[4], 3.0, 1e-15);
}
