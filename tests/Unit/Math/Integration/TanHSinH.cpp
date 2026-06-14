// SPDX-License-Identifier: Apache-2.0

#include "Math/Integration/TanHSinH.hpp"

#include <array>
#include <cmath>
#include <gtest/gtest.h>

TEST(MathIntegrationTanHSinH, IntegratesExponentialOnZeroToInfinity)
{
    const uv::math::integration::TanHSinH<double, 128> integrator;

    const double integral = integrator.integrateZeroToInf(
        [](double x)
        {
            return std::exp(-x);
        }
    );

    EXPECT_NEAR(integral, 1.0, 1e-12);
}

TEST(MathIntegrationTanHSinH, MultiIntegralMatchesScalarMoments)
{
    const uv::math::integration::TanHSinH<double, 128> integrator;

    const auto integrals = integrator.integrateZeroToInfMulti<2>(
        [](double x)
        {
            return std::array<double, 2>{std::exp(-x), x * std::exp(-x)};
        }
    );

    EXPECT_NEAR(integrals[0], 1.0, 1e-12);
    EXPECT_NEAR(integrals[1], 1.0, 1e-12);
}
