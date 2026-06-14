// SPDX-License-Identifier: Apache-2.0

#include "Core/Generate.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(IntegrationCoreGenerate, BuildsConsistentMarketState)
{
    const uv::core::MarketData<double> marketData{
        .interestRate = 0.04,
        .dividendYield = 0.01,
        .spot = 100.0
    };
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> moneyness{0.9, 1.0, 1.1};
    const uv::core::Matrix<double> vol{2, 3, 0.2};

    const auto state = uv::core::generateMarketState(
        marketData,
        std::span<const double>{maturities},
        std::span<const double>{moneyness},
        vol
    );

    EXPECT_NEAR(
        state.interestCurve.interpolateDF(1.0),
        std::exp(-marketData.interestRate),
        1e-15
    );
    EXPECT_NEAR(
        state.dividendCurve.interpolateDF(1.0),
        std::exp(-marketData.dividendYield),
        1e-15
    );
    EXPECT_NEAR(state.volSurface.forwards()[1], 100.0 * std::exp(0.03), 5e-14);
    EXPECT_DOUBLE_EQ(state.volSurface.strikes()[0], 90.0);
    EXPECT_DOUBLE_EQ(state.volSurface.vol()[1][2], 0.2);
}
