// SPDX-License-Identifier: Apache-2.0

#include "Golden.hpp"
#include "Math/Functions/Black.hpp"

#include <gtest/gtest.h>

TEST(RegressionMathBlack, AtmOneYearCallMatchesKnownValue)
{
    const auto golden =
        uv::tests::golden::readBlackKnownValue("tests/Golden/black_known_value.json");
    const double price = uv::math::black::priceB76(
        golden.t,
        golden.discountFactor,
        golden.forward,
        golden.volatility,
        golden.strike
    );

    EXPECT_NEAR(price, golden.price, golden.tolerance)
        << "t=" << golden.t << " discountFactor=" << golden.discountFactor
        << " forward=" << golden.forward << " volatility=" << golden.volatility
        << " strike=" << golden.strike << " actual=" << price
        << " expected=" << golden.price;
}
