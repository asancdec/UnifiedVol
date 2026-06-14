// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Math.hpp"

#include <gtest/gtest.h>

TEST(UnitModelsSVIMath, TotalVarianceMatchesRawSVIFormula)
{
    const double a = 0.02;
    const double b = 0.3;
    const double rho = -0.4;
    const double m = 0.1;
    const double sigma = 0.5;
    const double k = -0.2;

    const double x = k - m;
    const double expected = a + b * (rho * x + std::sqrt(x * x + sigma * sigma));

    EXPECT_NEAR(uv::models::svi::totalVariance(a, b, rho, m, sigma, k), expected, 1e-15);
}

TEST(UnitModelsSVIMath, AParameterMakesAtmVarianceMatch)
{
    const double atmVariance = 0.04;
    const double b = 0.3;
    const double rho = -0.4;
    const double m = 0.1;
    const double sigma = 0.5;

    const double a = uv::models::svi::detail::aParam(atmVariance, b, rho, m, sigma);

    EXPECT_NEAR(
        uv::models::svi::totalVariance(a, b, rho, m, sigma, 0.0),
        atmVariance,
        1e-15
    );
}

TEST(UnitModelsSVIMath, ButterflyDiagnosticIsPositiveForBenignParameters)
{
    EXPECT_GT(uv::models::svi::gk(0.04, 0.2, -0.3, 0.0, 0.5, 0.1), 0.0);
}
