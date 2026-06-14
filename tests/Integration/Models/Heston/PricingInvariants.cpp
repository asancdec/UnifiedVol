// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Price/Pricer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>

TEST(IntegrationHestonPricingInvariants, CallPricesRespectBoundsAndStrikeMonotonicity)
{
    uv::models::heston::price::Pricer<double, 160> pricer{};
    pricer.setParams({2.0, 0.05, 0.40, -0.60, 0.04});

    constexpr double t = 1.25;
    constexpr double dF = 0.97;
    constexpr double F = 100.0;
    const std::array<double, 7> strikes{60.0, 75.0, 90.0, 100.0, 110.0, 125.0, 140.0};

    double previous = dF * F;
    for (double K : strikes)
    {
        const double price = pricer.callPrice(t, dF, F, K);
        const double intrinsic = dF * std::max(F - K, 0.0);

        EXPECT_TRUE(std::isfinite(price)) << "K=" << K;
        EXPECT_GE(price + 5e-10, intrinsic) << "K=" << K;
        EXPECT_LE(price, dF * F + 5e-10) << "K=" << K;
        EXPECT_LE(price, previous + 5e-10) << "K=" << K;

        previous = price;
    }
}

TEST(IntegrationHestonPricingInvariants, CallPriceIncreasesWithForward)
{
    uv::models::heston::price::Pricer<double, 160> pricer{};
    pricer.setParams({1.8, 0.04, 0.35, -0.50, 0.04});

    constexpr double t = 1.0;
    constexpr double dF = 0.98;
    constexpr double K = 100.0;

    const double lowForward = pricer.callPrice(t, dF, 95.0, K);
    const double midForward = pricer.callPrice(t, dF, 100.0, K);
    const double highForward = pricer.callPrice(t, dF, 105.0, K);

    EXPECT_TRUE(std::isfinite(lowForward));
    EXPECT_TRUE(std::isfinite(midForward));
    EXPECT_TRUE(std::isfinite(highForward));
    EXPECT_LE(lowForward, midForward + 5e-10);
    EXPECT_LE(midForward, highForward + 5e-10);
}

TEST(IntegrationHestonPricingInvariants, SurfacePricesStayWithinNoArbitrageBounds)
{
    const std::array<double, 2> maturities{0.5, 1.5};
    const std::array<double, 2> discountFactors{0.99, 0.95};
    const std::array<double, 2> forwards{98.0, 103.0};
    const std::array<double, 5> strikes{70.0, 90.0, 100.0, 115.0, 140.0};

    uv::models::heston::price::Pricer<double, 160> pricer{};
    pricer.setParams({2.5, 0.04, 0.45, -0.70, 0.05});

    for (std::size_t i = 0; i < maturities.size(); ++i)
    {
        double previous = discountFactors[i] * forwards[i];
        for (double K : strikes)
        {
            const double price =
                pricer.callPrice(maturities[i], discountFactors[i], forwards[i], K);
            const double intrinsic = discountFactors[i] * std::max(forwards[i] - K, 0.0);

            EXPECT_TRUE(std::isfinite(price));
            EXPECT_GE(price + 5e-10, intrinsic);
            EXPECT_LE(price, discountFactors[i] * forwards[i] + 5e-10);
            EXPECT_LE(price, previous + 5e-10);

            previous = price;
        }
    }
}
