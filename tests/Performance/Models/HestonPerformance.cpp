// SPDX-License-Identifier: Apache-2.0

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Price/Pricer.hpp"
#include "Support/Performance/Budgets.hpp"
#include "Support/Performance/Timing.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(PerformanceHeston, PricesMediumSurfaceWithinThroughputBudget)
{
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        uv::tests::performance::HestonMediumSurfaceBudgetKey
    );
    const std::vector<double>
        maturities{0.08, 0.16, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0};
    const std::vector<double>
        forwards{100.0, 100.5, 101.0, 101.5, 102.0, 102.5, 103.0, 103.5, 104.0, 104.5};
    const std::vector<double> strikes{60.0,  65.0,  70.0,  75.0,  80.0,  85.0,  90.0,
                                      95.0,  100.0, 105.0, 110.0, 115.0, 120.0, 125.0,
                                      130.0, 135.0, 140.0, 145.0, 150.0, 160.0, 175.0};
    std::vector<double> moneyness;
    moneyness.reserve(strikes.size());
    for (const double strike : strikes)
        moneyness.emplace_back(strike / 100.0);

    const uv::core::Matrix<double> vols{maturities.size(), strikes.size(), 0.25};
    const uv::core::VolSurface<double>
        surface{maturities, forwards, strikes, moneyness, vols};
    const uv::core::Curve<double> curve{0.03, maturities};
    uv::models::heston::price::Pricer<double, 160> pricer{};
    pricer.setParams({2.0, 0.04, 0.35, -0.65, 0.05});

    auto prices = pricer.callPrice(surface, curve);
    ASSERT_EQ(prices.rows(), maturities.size());
    ASSERT_EQ(prices.cols(), strikes.size());

    const double ms = uv::tests::performance::bestElapsedMs(
        [&]
        {
            prices = pricer.callPrice(surface, curve);
        }
    );

    EXPECT_LT(ms, budget.maxMs);
    EXPECT_TRUE(std::isfinite(prices[0][0]));
    EXPECT_TRUE(std::isfinite(prices[prices.rows() - 1][prices.cols() - 1]));
}
