// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"

#include <cmath>
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

TEST(UnitModelsSVIMath, TotalVarianceParamsOverloadMatchesRawParameterOverload)
{
    const uv::models::svi::Params<double> params{1.0, 0.02, 0.3, -0.4, 0.1, 0.5};
    const double k = -0.2;

    EXPECT_NEAR(
        uv::models::svi::totalVariance(params, k),
        uv::models::svi::totalVariance(
            params.a,
            params.b,
            params.rho,
            params.m,
            params.sigma,
            k
        ),
        1e-15
    );
}

TEST(UnitModelsSVIMath, ButterflyDiagnosticIsPositiveForBenignParameters)
{
    EXPECT_GT(uv::models::svi::gk(0.04, 0.2, -0.3, 0.0, 0.5, 0.1), 0.0);
}

TEST(UnitModelsSVIMath, ButterflyDiagnosticParamsOverloadMatchesRawParameterOverload)
{
    const uv::models::svi::Params<double> params{1.0, 0.04, 0.2, -0.3, 0.0, 0.5};
    const double k = 0.1;

    EXPECT_DOUBLE_EQ(
        uv::models::svi::gk(params, k),
        uv::models::svi::gk(params.a, params.b, params.rho, params.m, params.sigma, k)
    );
}

TEST(UnitModelsSVIMath, ButterflyDiagnosticDetectsInvalidParameters)
{
    EXPECT_LT(uv::models::svi::gk(0.001, 5.0, 0.0, 0.0, 0.01, 0.5), 0.0);
}
