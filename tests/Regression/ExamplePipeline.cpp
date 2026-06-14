// SPDX-License-Identifier: Apache-2.0

#include "IO/Load.hpp"
#include "Models/Heston/BuildSurface.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

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

template <std::floating_point T> void expectSVIParamsNear(
    const uv::models::svi::Params<T>& actual,
    const uv::models::svi::Params<T>& expected,
    const T tolerance
)
{
    EXPECT_NEAR(actual.t, expected.t, tolerance);
    EXPECT_NEAR(actual.a, expected.a, tolerance);
    EXPECT_NEAR(actual.b, expected.b, tolerance);
    EXPECT_NEAR(actual.rho, expected.rho, tolerance);
    EXPECT_NEAR(actual.m, expected.m, tolerance);
    EXPECT_NEAR(actual.sigma, expected.sigma, tolerance);
}

template <std::floating_point T> void expectHestonParamsNear(
    const uv::models::heston::Params<T>& actual,
    const uv::models::heston::Params<T>& expected,
    const T tolerance
)
{
    EXPECT_NEAR(actual.kappa, expected.kappa, tolerance);
    EXPECT_NEAR(actual.theta, expected.theta, tolerance);
    EXPECT_NEAR(actual.sigma, expected.sigma, tolerance);
    EXPECT_NEAR(actual.rho, expected.rho, tolerance);
    EXPECT_NEAR(actual.v0, expected.v0, tolerance);
}
} // namespace

TEST(RegressionExamplePipeline, MainCppResultsRemainStable)
{
    using Real = double;
    constexpr Real tolerance = 1e-8;

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

    const uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP> sviOptimizer{
        uv::models::svi::detail::makeNLoptConfig(sviConfig)
    };
    const auto sviParams =
        uv::models::svi::calibrate(marketState.volSurface, sviOptimizer, false);
    const auto sviVolSurface =
        uv::models::svi::buildSurface(marketState.volSurface, sviParams);

    const auto hestonParams = uv::models::heston::calibrate::calibrate<Real>(
        sviVolSurface,
        marketState.interestCurve,
        hestonConfig
    );

    uv::models::heston::price::Pricer<Real> hestonPricer{};
    hestonPricer.setParams(hestonParams);
    const auto hestonVolSurface = uv::models::heston::buildSurface(
        sviVolSurface,
        marketState.interestCurve,
        hestonPricer
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

    const uv::Vector<uv::models::svi::Params<Real>> expectedSVIParams{
        {0.083333332999999996,
         -0.0454222144615542,
         0.149236855731602,
         0.07186134006880307,
         0.216970201378966,
         0.35679911238289619},
        {0.16666666699999999,
         -0.34981866687376706,
         1.1138461886956372,
         0.79558005431799084,
         0.93867382539223965,
         0.53374713287635056},
        {0.25,
         -0.32562925881393973,
         1.1176709797956048,
         0.78943538496968424,
         0.90125265042249558,
         0.49317518375943226},
        {0.33333333300000001,
         -0.3321410667384247,
         1.1241712488889342,
         0.77908837463747904,
         0.89702509048526358,
         0.49248616320897554},
        {0.41666666699999999,
         -0.34974817266854857,
         1.1311717936054102,
         0.76807803315653156,
         0.90421212190972311,
         0.50607232414172265},
        {0.5,
         -0.38037109273030467,
         1.1392419272638081,
         0.75555336591546474,
         0.92078057695565396,
         0.53420944189653552},
        {0.75,
         -0.33544940678689933,
         1.1381434904514793,
         0.7572476728805424,
         0.89705071489944921,
         0.48193506493710564},
        {1.0,
         -0.61994739213481764,
         1.1798007573948683,
         0.695201488441339,
         1.102542642233225,
         0.76660421138127122},
        {1.5,
         -0.46411964308195713,
         1.1574473349798731,
         0.72794039051001858,
         1.1141215410541001,
         0.63310633745924194},
        {2.0,
         -0.83331114115360672,
         1.1787248954911986,
         0.69674875592284291,
         1.6136221389666017,
         1.0351188825402371},
        {3.0,
         -8.8444781368502241,
         1.3080206169361941,
         -0.52902788694924741,
         -3.3612344508812764,
         7.9680550126940988}
    };
    const uv::models::heston::Params<Real> expectedHestonParams{
        11.166908774640188,
        0.061835473188871844,
        4.1454420263519616,
        -0.7480708082612596,
        0.37975758619313726
    };

    ASSERT_EQ(sviParams.size(), expectedSVIParams.size());
    for (std::size_t i = 0; i < sviParams.size(); ++i)
        expectSVIParamsNear(sviParams[i], expectedSVIParams[i], tolerance);
    expectHestonParamsNear(hestonParams, expectedHestonParams, tolerance);

    EXPECT_NEAR(marketState.volSurface.vol()[0][0], 1.0398700000000001, 1e-12);
    EXPECT_NEAR(marketState.volSurface.vol()[7][8], 0.24507000000000001, 1e-12);
    EXPECT_NEAR(marketState.volSurface.vol()[10][16], 0.15731999999999999, 1e-12);

    EXPECT_NEAR(sviVolSurface.vol()[0][0], 1.0432929310422308, tolerance);
    EXPECT_NEAR(sviVolSurface.vol()[7][8], 0.24506999999999993, tolerance);
    EXPECT_NEAR(sviVolSurface.vol()[10][16], 0.16051357088354287, tolerance);

    EXPECT_NEAR(hestonVolSurface.vol()[0][0], 0.97326102349139154, tolerance);
    EXPECT_NEAR(hestonVolSurface.vol()[7][8], 0.24270179257833455, tolerance);
    EXPECT_NEAR(hestonVolSurface.vol()[10][16], 0.1790726912635392, tolerance);

    const auto [meanError, maxError] =
        meanAndMaxAbsVolError(hestonVolSurface, sviVolSurface);

    EXPECT_NEAR(meanError, 0.011480071785533133, tolerance);
    EXPECT_NEAR(maxError, 0.070031907550839212, tolerance);
}
