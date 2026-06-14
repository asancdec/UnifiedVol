// SPDX-License-Identifier: Apache-2.0

#include "Math/Integration/TanHSinH.hpp"

#include <array>
#include <cmath>
#include <gtest/gtest.h>

TEST(IntegrationMathTanHSinH, IntegratesDampedMomentsOnPositiveHalfLine)
{
    const uv::math::integration::TanHSinH<double, 128> integrator;

    const auto moments = integrator.integrateZeroToInfMulti<3>(
        [](double x)
        {
            const double e = std::exp(-x);
            return std::array<double, 3>{e, x * e, x * x * e};
        }
    );

    EXPECT_NEAR(moments[0], 1.0, 1e-12);
    EXPECT_NEAR(moments[1], 1.0, 1e-12);
    EXPECT_NEAR(moments[2], 2.0, 1e-12);
}

TEST(IntegrationMathTanHSinH, IntegratesDampedOscillatoryFunction)
{
    const uv::math::integration::TanHSinH<double, 192> integrator;
    constexpr double a = 3.0;

    const double integral = integrator.integrateZeroToInf(
        [](double x)
        {
            return std::exp(-x) * std::cos(a * x);
        }
    );

    EXPECT_NEAR(integral, 1.0 / (1.0 + a * a), 1e-7);
}

TEST(IntegrationMathTanHSinH, MultiIntegralMatchesIndependentScalarRuns)
{
    const uv::math::integration::TanHSinH<double, 128> integrator;

    const auto multi = integrator.integrateZeroToInfMulti<3>(
        [](double x)
        {
            const double e = std::exp(-x);
            return std::array<double, 3>{
                e * std::cos(0.5 * x),
                e * std::cos(2.0 * x),
                e * std::cos(5.0 * x)
            };
        }
    );

    const std::array<double, 3> a{0.5, 2.0, 5.0};
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const double scalar = integrator.integrateZeroToInf(
            [ai = a[i]](double x)
            {
                return std::exp(-x) * std::cos(ai * x);
            }
        );
        EXPECT_NEAR(multi[i], scalar, 1e-14);
    }
}
