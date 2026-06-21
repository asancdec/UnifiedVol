// SPDX-License-Identifier: Apache-2.0

#include "IO/CSV/Load.hpp"
#include "Models/Heston/BuildSurface.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Optimization/NLopt/Optimizer.hpp"
#include "Support/Assertions.hpp"
#include "Support/Golden.hpp"

#include <cstddef>
#include <filesystem>
#include <gtest/gtest.h>

namespace
{
namespace assertions = uv::tests::assertions;
} // namespace

TEST(RegressionExamplePipeline, MainCppResultsRemainStable)
{
    using Real = double;

    const std::filesystem::path path{"data/VolSurface_SPY_04072011.csv"};
    const auto golden =
        uv::tests::golden::readExamplePipeline("tests/Golden/example_pipeline.json");
    const Real tolerance{golden.tolerance};
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

    const auto marketState = uv::io::csv::load::marketState(path, marketData);

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

    assertions::expectFiniteNonNegativeVolSurface(marketState.volSurface);
    assertions::expectFiniteNonNegativeVolSurface(sviVolSurface);
    assertions::expectFiniteNonNegativeVolSurface(hestonVolSurface);

    ASSERT_EQ(sviParams.size(), golden.sviParams.size());
    for (std::size_t i = 0; i < sviParams.size(); ++i)
        assertions::expectSVIParamsNear(sviParams[i], golden.sviParams[i], tolerance);
    assertions::expectHestonParamsNear(hestonParams, golden.hestonParams, tolerance);
    assertions::expectVolPointsNear(marketState.volSurface, golden.marketVols, tolerance);
    assertions::expectVolPointsNear(sviVolSurface, golden.sviVols, tolerance);
    assertions::expectVolPointsNear(hestonVolSurface, golden.hestonVols, tolerance);
    assertions::expectCallSurfaceNoArbitrage(
        hestonPricer.callPrice(sviVolSurface, marketState.interestCurve),
        sviVolSurface,
        marketState.interestCurve
    );

    const uv::tests::ErrorDiagnostics errors{
        assertions::volErrorDiagnostics(hestonVolSurface, sviVolSurface)
    };

    EXPECT_NEAR(errors.meanAbs, golden.meanAbsVolError, tolerance);
    EXPECT_NEAR(errors.maxAbs, golden.maxAbsVolError, tolerance);
    EXPECT_NEAR(errors.rmse, golden.rmseVolError, tolerance);
}
