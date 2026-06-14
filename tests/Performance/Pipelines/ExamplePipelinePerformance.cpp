// SPDX-License-Identifier: Apache-2.0

#include "IO/Load.hpp"
#include "Models/Heston/BuildSurface.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Support/Performance/Budgets.hpp"
#include "Support/Performance/Timing.hpp"

#include <filesystem>
#include <gtest/gtest.h>

TEST(PerformanceExamplePipeline, BuildsSVIAndHestonSurfacesWithinLatencyBudget)
{
    using Real = double;
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        uv::tests::performance::ExamplePipelineBudgetKey
    );

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
    std::size_t sviMaturities{};
    std::size_t hestonStrikes{};

    const double ms = uv::tests::performance::bestElapsedMs(
        [&]
        {
            const auto sviVolSurface =
                uv::models::svi::buildSurface(marketState, sviConfig);
            const auto hestonVolSurface = uv::models::heston::buildSurface<Real>(
                sviVolSurface,
                marketState.interestCurve,
                hestonConfig
            );
            sviMaturities = sviVolSurface.numMaturities();
            hestonStrikes = hestonVolSurface.numStrikes();
        }
    );

    EXPECT_EQ(sviMaturities, marketState.volSurface.numMaturities());
    EXPECT_EQ(hestonStrikes, marketState.volSurface.numStrikes());
    EXPECT_LT(ms, budget.maxMs);
}
