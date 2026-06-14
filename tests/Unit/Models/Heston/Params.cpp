// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Params.hpp"
#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(ModelsHestonParams, ConstructsFromSpanAndConvertsType)
{
    const std::vector<double> raw{1.0, 0.04, 0.3, -0.7, 0.05};

    const uv::models::heston::Params<long double> params{raw};
    const auto asDouble = params.as<double>();

    EXPECT_DOUBLE_EQ(asDouble.kappa, 1.0);
    EXPECT_DOUBLE_EQ(asDouble.theta, 0.04);
    EXPECT_DOUBLE_EQ(asDouble.sigma, 0.3);
    EXPECT_DOUBLE_EQ(asDouble.rho, -0.7);
    EXPECT_DOUBLE_EQ(asDouble.v0, 0.05);
}

TEST(ModelsHestonParams, RejectsWrongSpanSize)
{
    const std::vector<double> raw{1.0, 0.04, 0.3};

    EXPECT_THROW((uv::models::heston::Params<double>{raw}), uv::errors::UnifiedVolError);
}
