// SPDX-License-Identifier: Apache-2.0

#include "Budgets.hpp"
#include "Core/Matrix.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"
#include "Timing.hpp"

#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(PerformanceSVI, CalibratesSyntheticSurfaceWithinLatencyBudget)
{
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        "sviSyntheticCalibration"
    );
    const std::vector<double> maturities{0.08, 0.25, 0.5, 1.0, 2.0};
    const std::vector<double> logKF{-0.45, -0.30, -0.15, 0.0, 0.15, 0.30, 0.45};
    const uv::Vector<uv::models::svi::Params<double>> truth{
        {0.08, 0.010, 0.22, -0.70, -0.03, 0.22},
        {0.25, 0.025, 0.26, -0.65, -0.02, 0.26},
        {0.5, 0.045, 0.30, -0.55, -0.01, 0.30},
        {1.0, 0.080, 0.34, -0.45, 0.00, 0.36},
        {2.0, 0.140, 0.38, -0.35, 0.01, 0.44}
    };

    uv::core::Matrix<double> logKFMatrix{maturities.size(), logKF.size()};
    uv::core::Matrix<double> totalVariance{maturities.size(), logKF.size()};
    for (std::size_t i = 0; i < maturities.size(); ++i)
    {
        for (std::size_t j = 0; j < logKF.size(); ++j)
        {
            logKFMatrix[i][j] = logKF[j];
            const auto& p = truth[i];
            totalVariance[i][j] =
                uv::models::svi::totalVariance(p.a, p.b, p.rho, p.m, p.sigma, logKF[j]);
        }
    }

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
                std::span<const double>{maturities},
                logKFMatrix,
                totalVariance,
                optimizer
            );
        }
    );

    EXPECT_EQ(calibrated.size(), maturities.size());
    EXPECT_LT(ms, budget.maxMs);
}
