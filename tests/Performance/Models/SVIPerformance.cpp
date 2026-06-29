// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/Config.hpp"
#include "Optimization/NLopt/Optimizer.hpp"
#include "Support/Models/SVI.hpp"
#include "Support/Performance/Budgets.hpp"
#include "Support/Performance/Timing.hpp"

#include <gtest/gtest.h>
#include <span>

TEST(PerformanceSVI, CalibratesSyntheticSurfaceWithinLatencyBudget)
{
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        uv::tests::performance::SVISyntheticCalibrationBudgetKey
    );
    const auto data = uv::tests::models::svi::makePerformanceSyntheticSurfaceCase();

    const uv::models::svi::Config config{
        .objectiveTol = 1e-12,
        .maxEval = 2000,
        .verbose = false,
        .printParams = false
    };
    const uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        uv::models::svi::detail::makeNLoptConfig(config)
    };
    uv::Vector<uv::models::svi::Params<double>> calibrated;

    const double ms = uv::tests::performance::bestElapsedMs(
        [&]
        {
            calibrated = uv::models::svi::calibrate(
                std::span<const double>{data.maturities},
                data.logKFMatrix,
                data.totalVariance,
                optimizer
            );
        }
    );

    EXPECT_EQ(calibrated.size(), data.maturities.size());
    EXPECT_LT(ms, budget.maxMs);
}
