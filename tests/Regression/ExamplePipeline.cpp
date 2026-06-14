// SPDX-License-Identifier: Apache-2.0

#include "IO/Load.hpp"
#include "Models/Heston/BuildSurface.hpp"
#include "Models/SVI/BuildSurface.hpp"

#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>

namespace
{
template <std::floating_point T>
void expectFiniteNonNegativeVolSurface(const uv::core::VolSurface<T>& surface)
{
    for (std::size_t i = 0; i < surface.numMaturities(); ++i)
    {
        for (std::size_t j = 0; j < surface.numStrikes(); ++j)
        {
            EXPECT_TRUE(std::isfinite(surface.vol()[i][j]));
            EXPECT_GE(surface.vol()[i][j], T{0});
        }
    }
}

template <std::floating_point T> std::pair<T, T> meanAndMaxAbsVolError(
    const uv::core::VolSurface<T>& lhs,
    const uv::core::VolSurface<T>& rhs
)
{
    T sum{0};
    T maxError{0};
    std::size_t count{0};

    for (std::size_t i = 0; i < lhs.numMaturities(); ++i)
    {
        for (std::size_t j = 0; j < lhs.numStrikes(); ++j)
        {
            const T error = std::abs(lhs.vol()[i][j] - rhs.vol()[i][j]);
            sum += error;
            maxError = std::max(maxError, error);
            ++count;
        }
    }

    return {sum / static_cast<T>(count), maxError};
}
} // namespace

TEST(RegressionExamplePipeline, MainCppResultsRemainStable)
{
    using Real = double;

    const std::filesystem::path path{"data/VolSurface_SPY_04072011.csv"};
    const uv::core::MarketData<Real> marketData{
        .interestRate = 0.0,
        .dividendYield = 0.0,
        .spot = 485.77548
    };
    const uv::models::svi::Config sviConfig{
        .objectiveTol = 1e-12,
        .maxEval = 10000,
        .verbose = false,
        .printParams = false
    };
    const uv::models::heston::calibrate::Config hestonConfig{
        .tolerance = 1e-11,
        .maxEval = 10000,
        .verbosity = uv::opt::ceres::Verbosity::None
    };

    const auto marketState = uv::io::load::marketState(path, marketData);
    const auto sviVolSurface = uv::models::svi::buildSurface(marketState, sviConfig);
    const auto hestonVolSurface = uv::models::heston::buildSurface<Real>(
        sviVolSurface,
        marketState.interestCurve,
        hestonConfig
    );

    ASSERT_EQ(marketState.volSurface.numMaturities(), 11U);
    ASSERT_EQ(marketState.volSurface.numStrikes(), 17U);
    ASSERT_EQ(sviVolSurface.numMaturities(), marketState.volSurface.numMaturities());
    ASSERT_EQ(sviVolSurface.numStrikes(), marketState.volSurface.numStrikes());
    ASSERT_EQ(hestonVolSurface.numMaturities(), marketState.volSurface.numMaturities());
    ASSERT_EQ(hestonVolSurface.numStrikes(), marketState.volSurface.numStrikes());

    expectFiniteNonNegativeVolSurface(marketState.volSurface);
    expectFiniteNonNegativeVolSurface(sviVolSurface);
    expectFiniteNonNegativeVolSurface(hestonVolSurface);

    EXPECT_NEAR(marketState.volSurface.vol()[0][0], 1.03987, 1e-12);
    EXPECT_NEAR(marketState.volSurface.vol()[7][8], 0.24507, 1e-12);
    EXPECT_NEAR(marketState.volSurface.vol()[10][16], 0.15732, 1e-12);

    EXPECT_NEAR(sviVolSurface.vol()[0][0], 1.04329, 5e-5);
    EXPECT_NEAR(sviVolSurface.vol()[7][8], 0.24507, 5e-5);
    EXPECT_NEAR(sviVolSurface.vol()[10][16], 0.16051, 5e-5);

    EXPECT_NEAR(hestonVolSurface.vol()[0][0], 0.97326, 5e-5);
    EXPECT_NEAR(hestonVolSurface.vol()[7][8], 0.24270, 5e-5);
    EXPECT_NEAR(hestonVolSurface.vol()[10][16], 0.17907, 5e-5);

    const auto [meanError, maxError] =
        meanAndMaxAbsVolError(hestonVolSurface, sviVolSurface);

    EXPECT_LT(meanError, 0.025);
    EXPECT_LT(maxError, 0.125);
}
