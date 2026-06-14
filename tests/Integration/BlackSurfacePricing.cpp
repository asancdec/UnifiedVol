// SPDX-License-Identifier: Apache-2.0

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Functions/Black.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(IntegrationMathBlack, PricesWholeVolSurfaceWithCurve)
{
    const std::vector<double> maturities{1.0, 2.0};
    const std::vector<double> forwards{100.0, 102.0};
    const std::vector<double> strikes{90.0, 100.0};
    const std::vector<double> moneyness{0.9, 1.0};
    uv::core::Matrix<double> vols{2, 2};
    vols[0][0] = 0.20;
    vols[0][1] = 0.21;
    vols[1][0] = 0.22;
    vols[1][1] = 0.23;

    const uv::core::VolSurface<double>
        surface{maturities, forwards, strikes, moneyness, vols};
    const uv::core::Curve<double> curve{0.03, maturities};

    const auto prices = uv::math::black::priceB76(surface, curve);

    EXPECT_EQ(prices.rows(), surface.numMaturities());
    EXPECT_EQ(prices.cols(), surface.numStrikes());
    EXPECT_DOUBLE_EQ(
        prices[0][1],
        uv::math::black::priceB76(1.0, curve.interpolateDF(1.0), 100.0, 0.21, 100.0)
    );
    EXPECT_DOUBLE_EQ(
        prices[1][0],
        uv::math::black::priceB76(2.0, curve.interpolateDF(2.0), 102.0, 0.22, 90.0)
    );
}
