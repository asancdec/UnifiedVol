// SPDX-License-Identifier: Apache-2.0

#include "Core/VolSurface.hpp"
#include "Base/Errors/Errors.hpp"
#include "Core/Matrix.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(CoreVolSurface, StoresGridAndVolMatrix)
{
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> forwards{100.0, 101.0};
    const std::vector<double> strikes{90.0, 100.0, 110.0};
    const std::vector<double> moneyness{0.9, 1.0, 1.1};
    uv::core::Matrix<double> vols{2, 3};
    vols[0][0] = 0.20;
    vols[0][1] = 0.21;
    vols[0][2] = 0.22;
    vols[1][0] = 0.23;
    vols[1][1] = 0.24;
    vols[1][2] = 0.25;

    const uv::core::VolSurface<double>
        surface{maturities, forwards, strikes, moneyness, vols};

    EXPECT_EQ(surface.numMaturities(), 2U);
    EXPECT_EQ(surface.numStrikes(), 3U);
    EXPECT_DOUBLE_EQ(surface.maturities()[1], 1.0);
    EXPECT_DOUBLE_EQ(surface.forwards()[0], 100.0);
    EXPECT_DOUBLE_EQ(surface.strikes()[2], 110.0);
    EXPECT_DOUBLE_EQ(surface.moneyness()[0], 0.9);
    EXPECT_DOUBLE_EQ(surface.vol()[1][2], 0.25);
}

TEST(CoreVolSurface, RejectsShapeMismatch)
{
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> forwards{100.0, 101.0};
    const std::vector<double> strikes{90.0, 100.0};
    const std::vector<double> moneyness{0.9, 1.0};
    const uv::core::Matrix<double> wrongShape{1, 2, 0.2};

    EXPECT_THROW(
        (uv::core::VolSurface<
            double>{maturities, forwards, strikes, moneyness, wrongShape}),
        uv::errors::UnifiedVolError
    );
}

TEST(CoreVolSurface, RejectsZeroForwardStrikeAndMoneyness)
{
    const std::vector<double> maturities{0.5, 1.0};
    const std::vector<double> forwards{100.0, 101.0};
    const std::vector<double> strikes{90.0, 100.0};
    const std::vector<double> moneyness{0.9, 1.0};
    const uv::core::Matrix<double> vols{2, 2, 0.2};

    {
        auto badForwards{forwards};
        badForwards[0] = 0.0;
        EXPECT_THROW(
            (uv::core::VolSurface<
                double>{maturities, badForwards, strikes, moneyness, vols}),
            uv::errors::UnifiedVolError
        );
    }

    {
        auto badStrikes{strikes};
        badStrikes[0] = 0.0;
        EXPECT_THROW(
            (uv::core::VolSurface<
                double>{maturities, forwards, badStrikes, moneyness, vols}),
            uv::errors::UnifiedVolError
        );
    }

    {
        auto badMoneyness{moneyness};
        badMoneyness[0] = 0.0;
        EXPECT_THROW(
            (uv::core::VolSurface<
                double>{maturities, forwards, strikes, badMoneyness, vols}),
            uv::errors::UnifiedVolError
        );
    }
}
