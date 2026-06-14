// SPDX-License-Identifier: Apache-2.0

#include "IO/Load.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"

#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <vector>

TEST(IntegrationLoadAndReport, LoadsExampleCsv)
{
    const std::filesystem::path path{"data/VolSurface_SPY_04072011.csv"};
    const uv::core::MarketData<double> marketData{
        .interestRate = 0.01,
        .dividendYield = 0.02,
        .spot = 485.77548
    };

    const auto marketState = uv::io::load::marketState(path, marketData);

    ASSERT_EQ(marketState.volSurface.numMaturities(), 11U);
    ASSERT_EQ(marketState.volSurface.numStrikes(), 17U);
    EXPECT_DOUBLE_EQ(marketState.volSurface.maturities()[0], 0.083333333);
    EXPECT_DOUBLE_EQ(marketState.volSurface.moneyness()[16], 1.5);
    EXPECT_DOUBLE_EQ(marketState.volSurface.vol()[0][0], 1.03987);
}

TEST(IntegrationLoadAndDerivedResults, ComputesDerivedResultsForMarketState)
{
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> forwards{100.0, 101.0};
    const std::vector<double> strikes{90.0, 100.0};
    const std::vector<double> moneyness{0.9, 1.0};
    const uv::core::Matrix<double> vols{2, 2, 0.2};
    const uv::core::MarketState<double> marketState{
        .interestCurve = uv::core::Curve<double>{0.03, maturities},
        .dividendCurve = uv::core::Curve<double>{0.01, maturities},
        .volSurface =
            uv::core::VolSurface<double>{maturities, forwards, strikes, moneyness, vols}
    };

    const auto variance = uv::math::vol::variance(marketState.volSurface);
    const auto totalVariance = uv::math::vol::totalVariance(marketState.volSurface);
    const auto logKF = uv::math::vol::logKF(marketState.volSurface);
    const auto calls = uv::math::black::priceB76(
        marketState.volSurface,
        marketState.interestCurve,
        true
    );
    const auto puts = uv::math::black::priceB76(
        marketState.volSurface,
        marketState.interestCurve,
        false
    );

    EXPECT_DOUBLE_EQ(variance[0][0], 0.04);
    EXPECT_DOUBLE_EQ(totalVariance[0][0], 0.02);
    EXPECT_DOUBLE_EQ(totalVariance[1][1], 0.04);
    EXPECT_NEAR(logKF[0][0], std::log(90.0 / 100.0), 1e-15);
    EXPECT_NEAR(logKF[1][1], std::log(100.0 / 101.0), 1e-15);

    EXPECT_GT(calls[0][0], calls[0][1]);
    EXPECT_LT(puts[0][0], puts[0][1]);
    EXPECT_NEAR(
        calls[1][1] - puts[1][1],
        marketState.interestCurve.interpolateDF(1.0) * (101.0 - 100.0),
        1e-12
    );
}
