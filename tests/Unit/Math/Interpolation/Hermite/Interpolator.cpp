// SPDX-License-Identifier: Apache-2.0

#include "Math/Interpolation/Hermite/Interpolator.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(MathHermiteInterpolator, ReproducesLinearData)
{
    const std::vector<double> xs{0.0, 1.0, 2.0, 3.0};
    const std::vector<double> ys{1.0, 3.0, 5.0, 7.0};
    const std::vector<double> x{0.5, 1.5, 2.5};

    const auto y = uv::math::interp::hermite::PchipInterpolator<double>{}(x, xs, ys);

    ASSERT_EQ(y.size(), x.size());
    EXPECT_NEAR(y[0], 2.0, 1e-15);
    EXPECT_NEAR(y[1], 4.0, 1e-15);
    EXPECT_NEAR(y[2], 6.0, 1e-15);
}

TEST(MathHermiteInterpolator, UsesProvidedDerivatives)
{
    const std::vector<double> xs{0.0, 1.0};
    const std::vector<double> ys{0.0, 1.0};
    const std::vector<double> dydx{0.0, 0.0};
    const std::vector<double> x{0.25, 0.5, 0.75};

    const auto y =
        uv::math::interp::hermite::PchipInterpolator<double>{}(x, xs, ys, dydx);

    EXPECT_NEAR(y[0], 0.15625, 1e-15);
    EXPECT_NEAR(y[1], 0.5, 1e-15);
    EXPECT_NEAR(y[2], 0.84375, 1e-15);
}

TEST(MathHermiteInterpolator, ClampsOutsideInputDomain)
{
    const std::vector<double> xs{1.0, 2.0, 3.0};
    const std::vector<double> ys{10.0, 20.0, 30.0};

    const auto yLeft =
        uv::math::interp::hermite::PchipInterpolator<double>{}(0.0, xs, ys);
    const auto yRight =
        uv::math::interp::hermite::PchipInterpolator<double>{}(4.0, xs, ys);

    EXPECT_DOUBLE_EQ(yLeft, 10.0);
    EXPECT_DOUBLE_EQ(yRight, 30.0);
}
