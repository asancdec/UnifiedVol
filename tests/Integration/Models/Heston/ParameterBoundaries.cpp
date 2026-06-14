// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <array>
#include <cmath>
#include <gtest/gtest.h>

TEST(IntegrationHestonParameterBoundaries, NearBoundaryParametersProduceFinitePrices)
{
    uv::models::heston::price::Pricer<double, 192> pricer{};
    const std::array<uv::models::heston::Params<double>, 3> params{
        uv::models::heston::Params<double>{8.0, 0.02, 1e-8, -0.999, 0.02},
        uv::models::heston::Params<double>{1.2, 0.20, 1.00, 0.999, 0.25},
        uv::models::heston::Params<double>{0.05, 0.01, 0.80, -0.95, 0.30}
    };

    for (const auto& p : params)
    {
        pricer.setParams(p);
        const double price = pricer.callPrice(1.0, 0.98, 100.0, 100.0);
        EXPECT_TRUE(std::isfinite(price));
        EXPECT_GE(price, -1e-8);
        EXPECT_LE(price, 98.0 + 1e-8);
    }
}

TEST(IntegrationHestonPricingValidation, RejectsInvalidCallInputs)
{
    uv::models::heston::price::Pricer<double, 64> pricer{};
    pricer.setParams({2.0, 0.04, 0.35, -0.6, 0.04});

    EXPECT_THROW(pricer.callPrice(1.0, 0.0, 100.0, 100.0), uv::errors::UnifiedVolError);
    EXPECT_THROW(pricer.callPrice(1.0, 0.98, -100.0, 100.0), uv::errors::UnifiedVolError);
    EXPECT_THROW(pricer.callPrice(1.0, 0.98, 100.0, -100.0), uv::errors::UnifiedVolError);
}
