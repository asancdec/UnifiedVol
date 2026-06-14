// SPDX-License-Identifier: Apache-2.0

#include "Math/Functions/Black.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(MathBlack, B76CallPutParityHolds)
{
    const double t = 1.5;
    const double dF = 0.96;
    const double F = 105.0;
    const double K = 100.0;
    const double vol = 0.22;

    const double call = uv::math::black::priceB76(t, dF, F, vol, K);
    const double put = uv::math::black::priceB76(t, dF, F, vol, K, true, false);

    EXPECT_NEAR(call - put, dF * (F - K), 1e-13);
}

TEST(MathBlack, BlackScholesMatchesBlack76ForwardForm)
{
    const double t = 2.0;
    const double r = 0.04;
    const double q = 0.01;
    const double S = 100.0;
    const double K = 110.0;
    const double vol = 0.25;
    const double dF = std::exp(-r * t);
    const double F = S * std::exp((r - q) * t);

    const double bs = uv::math::black::priceBS(t, r, q, vol, S, K);
    const double b76 = uv::math::black::priceB76(t, dF, F, vol, K);

    EXPECT_NEAR(bs, b76, 1e-12);
}

TEST(MathBlack, VectorizedB76MatchesScalarPrices)
{
    const std::vector<double> vols{0.20, 0.25, 0.30};
    const std::vector<double> strikes{90.0, 100.0, 110.0};
    std::vector<double> out(3);

    uv::math::black::priceB76(
        std::span<double>{out},
        1.0,
        0.98,
        100.0,
        std::span<const double>{vols},
        std::span<const double>{strikes}
    );

    for (std::size_t i = 0; i < out.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(
            out[i],
            uv::math::black::priceB76(1.0, 0.98, 100.0, vols[i], strikes[i])
        );
    }
}

TEST(MathBlack, VegaIsPositiveForRegularOption)
{
    EXPECT_GT(uv::math::black::vegaB76(1.0, 0.98, 100.0, 0.2, 100.0), 0.0);
}
