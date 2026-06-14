// SPDX-License-Identifier: Apache-2.0

#include "../../../Support/Tolerances.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>

namespace
{
constexpr double StressNoArbTolerance = 2e-6;

struct HestonStressScenario
{
    uv::models::heston::Params<double> params;
    double t{};
    double discountFactor{};
    double forward{};
};

void expectFiniteBoundedMonotoneCalls(
    uv::models::heston::price::Pricer<double, 192>& pricer,
    const HestonStressScenario& scenario,
    const std::array<double, 7>& strikes
)
{
    pricer.setParams(scenario.params);

    double previous = scenario.discountFactor * scenario.forward;
    for (const double strike : strikes)
    {
        const double price =
            pricer
                .callPrice(scenario.t, scenario.discountFactor, scenario.forward, strike);
        const double intrinsic =
            scenario.discountFactor * std::max(scenario.forward - strike, 0.0);

        EXPECT_TRUE(std::isfinite(price)) << "t=" << scenario.t << " K=" << strike;
        EXPECT_GE(price + StressNoArbTolerance, intrinsic)
            << "t=" << scenario.t << " K=" << strike;
        EXPECT_LE(
            price,
            scenario.discountFactor * scenario.forward + StressNoArbTolerance
        ) << "t="
          << scenario.t << " K=" << strike;
        EXPECT_LE(price, previous + StressNoArbTolerance)
            << "t=" << scenario.t << " K=" << strike;
        previous = price;
    }
}
} // namespace

TEST(StressHestonPricing, HandlesExtremeMaturitiesAndCorrelations)
{
    uv::models::heston::price::Pricer<double, 192> pricer{};
    const std::array<double, 7> strikes{40.0, 60.0, 80.0, 100.0, 120.0, 160.0, 220.0};
    const std::array<HestonStressScenario, 4> scenarios{
        HestonStressScenario{{4.0, 0.02, 0.20, -0.95, 0.01}, 1.0 / 365.0, 0.999, 100.0},
        HestonStressScenario{{1.0, 0.20, 1.20, -0.90, 0.25}, 0.05, 0.998, 100.0},
        HestonStressScenario{{0.35, 0.09, 0.80, 0.90, 0.16}, 10.0, 0.75, 100.0},
        HestonStressScenario{{8.0, 0.04, 0.10, -0.10, 0.25}, 5.0, 0.85, 125.0}
    };

    for (const auto& scenario : scenarios)
        expectFiniteBoundedMonotoneCalls(pricer, scenario, strikes);
}
