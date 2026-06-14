// SPDX-License-Identifier: Apache-2.0

#include "Math/Functions/Black.hpp"

#include <gtest/gtest.h>

TEST(RegressionMathBlack, AtmOneYearCallMatchesKnownValue)
{
    const double price = uv::math::black::priceB76(1.0, 0.95, 100.0, 0.20, 100.0);

    EXPECT_NEAR(price, 7.567289082635459, 1e-12);
}
