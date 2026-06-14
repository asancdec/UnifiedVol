// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Params.hpp"
#include "Models/SVI/Math.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(ModelsSVIParams, StoresRawParamsAndConvertsType)
{
    const uv::models::svi::Params<double> params{1.0, 0.02, 0.3, -0.4, 0.1, 0.5};
    const auto asFloat = params.as<float>();

    EXPECT_FLOAT_EQ(asFloat.t, 1.0F);
    EXPECT_FLOAT_EQ(asFloat.a, 0.02F);
    EXPECT_FLOAT_EQ(asFloat.b, 0.3F);
    EXPECT_FLOAT_EQ(asFloat.rho, -0.4F);
    EXPECT_FLOAT_EQ(asFloat.m, 0.1F);
    EXPECT_FLOAT_EQ(asFloat.sigma, 0.5F);
}

TEST(ModelsSVIParams, ComputesAFromAtmTotalVarianceConstructor)
{
    const std::vector<double> packed{0.3, -0.4, 0.1, 0.5};
    const double atmTotalVariance = 0.04;

    const uv::models::svi::Params<double> params{1.0, packed, atmTotalVariance};

    EXPECT_DOUBLE_EQ(params.t, 1.0);
    EXPECT_DOUBLE_EQ(params.b, packed[0]);
    EXPECT_DOUBLE_EQ(params.rho, packed[1]);
    EXPECT_DOUBLE_EQ(params.m, packed[2]);
    EXPECT_DOUBLE_EQ(params.sigma, packed[3]);
    EXPECT_NEAR(
        uv::models::svi::totalVariance(
            params.a,
            params.b,
            params.rho,
            params.m,
            params.sigma,
            0.0
        ),
        atmTotalVariance,
        1e-15
    );
}
