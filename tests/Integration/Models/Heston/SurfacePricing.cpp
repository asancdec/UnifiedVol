// SPDX-License-Identifier: Apache-2.0

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(IntegrationHestonSurfacePricing, PricesWholeSurfaceWithStoredParams)
{
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> forwards{100.0, 102.0};
    const std::vector<double> strikes{90.0, 100.0, 110.0};
    const std::vector<double> moneyness{0.9, 1.0, 1.1};
    const uv::core::Matrix<double> vols{2, 3, 0.2};
    const uv::core::VolSurface<double>
        surface{maturities, forwards, strikes, moneyness, vols};
    const uv::core::Curve<double> curve{0.03, maturities};
    uv::models::heston::price::Pricer<double, 64> pricer{};
    pricer.setParams({2.0, 0.04, 0.35, -0.6, 0.04});

    const auto prices = pricer.callPrice(surface, curve);

    ASSERT_EQ(prices.rows(), surface.numMaturities());
    ASSERT_EQ(prices.cols(), surface.numStrikes());
    for (std::size_t i = 0; i < prices.rows(); ++i)
    {
        for (std::size_t j = 0; j < prices.cols(); ++j)
        {
            EXPECT_TRUE(std::isfinite(prices[i][j]));
            EXPECT_GE(prices[i][j], 0.0);
        }
    }
}

TEST(IntegrationHestonSurfacePricing, StoredParamsValidatePricingInputs)
{
    uv::models::heston::price::Pricer<double, 64> pricer{};

    EXPECT_THROW(pricer.callPrice(1.0, 0.98, 100.0, 100.0), uv::errors::UnifiedVolError);

    pricer.setParams({2.0, 0.04, 0.35, -0.6, 0.04});

    EXPECT_THROW(pricer.callPrice(1.0, 0.98, 0.0, 100.0), uv::errors::UnifiedVolError);
    EXPECT_THROW(pricer.callPrice(1.0, 0.98, 100.0, 0.0), uv::errors::UnifiedVolError);
    EXPECT_THROW(pricer.callPrice(0.0, 0.98, 100.0, 100.0), uv::errors::UnifiedVolError);
}

TEST(IntegrationHestonSurfacePricing, RejectsInvalidAlphaConfig)
{
    auto quad = std::make_shared<const uv::math::integration::TanHSinH<double, 64>>();

    EXPECT_THROW(
        (uv::models::heston::price::Pricer<double, 64>{
            quad,
            {.alphaItm = 1.0, .alphaOtm = 2.0}
        }),
        uv::errors::UnifiedVolError
    );

    EXPECT_THROW(
        (uv::models::heston::price::Pricer<double, 64>{
            quad,
            {.alphaItm = -2.0, .alphaOtm = 0.0}
        }),
        uv::errors::UnifiedVolError
    );
}
