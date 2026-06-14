// SPDX-License-Identifier: Apache-2.0

#include "Math/Functions/Volatility.hpp"
#include "Math/Functions/Black.hpp"

#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"

#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(MathVolatility, LogKfVarianceAndTotalVarianceAreComputedElementwise)
{
    const std::vector<double> strikes{90.0, 100.0, 110.0};
    std::vector<double> logKF(3);
    std::vector<double> variances(3);
    std::vector<double> totalVariances(3);
    const std::vector<double> vols{0.20, 0.25, 0.30};

    uv::math::vol::logKF(
        std::span<double>{logKF},
        100.0,
        std::span<const double>{strikes}
    );
    uv::math::vol::variance(std::span<double>{variances}, std::span<const double>{vols});
    uv::math::vol::totalVariance(
        std::span<double>{totalVariances},
        2.0,
        std::span<const double>{vols}
    );

    EXPECT_NEAR(logKF[0], std::log(0.9), 1e-15);
    EXPECT_DOUBLE_EQ(logKF[1], 0.0);
    EXPECT_DOUBLE_EQ(variances[0], 0.04);
    EXPECT_DOUBLE_EQ(variances[2], 0.09);
    EXPECT_DOUBLE_EQ(totalVariances[0], 0.08);
    EXPECT_DOUBLE_EQ(totalVariances[2], 0.18);
}

TEST(MathVolatility, VolFromTotalVarianceInvertsTotalVariance)
{
    const std::vector<double> vols{0.20, 0.35};
    std::vector<double> totalVariances(2);
    std::vector<double> recovered(2);

    uv::math::vol::totalVariance(
        std::span<double>{totalVariances},
        3.0,
        std::span<const double>{vols}
    );
    uv::math::vol::volFromTotalVariance(
        std::span<double>{recovered},
        3.0,
        std::span<const double>{totalVariances}
    );

    EXPECT_NEAR(recovered[0], vols[0], 1e-15);
    EXPECT_NEAR(recovered[1], vols[1], 1e-15);
}

TEST(MathVolatility, ImpliedVolRoundTripsBlack76CallPrice)
{
    const double t = 1.25;
    const double dF = 0.97;
    const double F = 100.0;
    const double K = 105.0;
    const double vol = 0.23;
    const double call = uv::math::black::priceB76(t, dF, F, vol, K);

    const double implied = uv::math::vol::impliedVol(call, t, dF, F, K);

    EXPECT_NEAR(implied, vol, 1e-12);
}

TEST(MathVolatility, SurfaceOverloadsComputeLogKfVarianceAndTotalVariance)
{
    const std::vector<double> maturities{1.0, 2.0};
    const std::vector<double> forwards{100.0, 110.0};
    const std::vector<double> strikes{90.0, 100.0};
    const std::vector<double> moneyness{0.9, 1.0};
    uv::core::Matrix<double> vols{2, 2};
    vols[0][0] = 0.20;
    vols[0][1] = 0.25;
    vols[1][0] = 0.30;
    vols[1][1] = 0.35;
    const uv::core::VolSurface<double>
        surface{maturities, forwards, strikes, moneyness, vols};

    const auto logKF = uv::math::vol::logKF(surface);
    const auto variance = uv::math::vol::variance(surface);
    const auto totalVariance = uv::math::vol::totalVariance(surface);

    EXPECT_NEAR(logKF[0][0], std::log(0.9), 1e-15);
    EXPECT_NEAR(logKF[1][1], std::log(100.0 / 110.0), 1e-15);
    EXPECT_DOUBLE_EQ(variance[0][1], 0.0625);
    EXPECT_DOUBLE_EQ(totalVariance[1][0], 0.18);
}

TEST(MathVolatility, MatrixVolFromTotalVarianceInvertsSurfaceTotalVariance)
{
    const std::vector<double> maturities{1.0, 2.0};
    uv::core::Matrix<double> totalVariance{2, 2};
    totalVariance[0][0] = 0.04;
    totalVariance[0][1] = 0.09;
    totalVariance[1][0] = 0.18;
    totalVariance[1][1] = 0.32;

    const auto vol = uv::math::vol::volFromTotalVariance(
        std::span<const double>{maturities},
        totalVariance
    );

    EXPECT_NEAR(vol[0][0], 0.20, 1e-15);
    EXPECT_NEAR(vol[0][1], 0.30, 1e-15);
    EXPECT_NEAR(vol[1][0], 0.30, 1e-15);
    EXPECT_NEAR(vol[1][1], 0.40, 1e-15);
}

TEST(MathVolatility, AtmParameterInterpolatesAtZeroLogMoneyness)
{
    const std::vector<double> logKF{-0.2, 0.0, 0.2};
    const std::vector<double> parameters{0.5, 1.0, 1.5};

    EXPECT_NEAR(
        uv::math::vol::atmParameter(
            std::span<const double>{parameters},
            std::span<const double>{logKF}
        ),
        1.0,
        1e-15
    );
}

TEST(MathVolatility, MatrixImpliedVolRoundTripsBlackSurfacePrices)
{
    const std::vector<double> maturities{1.0, 2.0};
    const std::vector<double> discountFactors{0.98, 0.95};
    const std::vector<double> forwards{100.0, 105.0};
    const std::vector<double> strikes{90.0, 100.0};
    uv::core::Matrix<double> callPrices{2, 2};

    for (std::size_t i = 0; i < maturities.size(); ++i)
    {
        for (std::size_t j = 0; j < strikes.size(); ++j)
        {
            callPrices[i][j] = uv::math::black::priceB76(
                maturities[i],
                discountFactors[i],
                forwards[i],
                0.25,
                strikes[j]
            );
        }
    }

    const auto implied = uv::math::vol::impliedVol(
        callPrices,
        std::span<const double>{maturities},
        std::span<const double>{discountFactors},
        std::span<const double>{forwards},
        std::span<const double>{strikes}
    );

    EXPECT_NEAR(implied[0][0], 0.25, 1e-12);
    EXPECT_NEAR(implied[0][1], 0.25, 1e-12);
    EXPECT_NEAR(implied[1][0], 0.25, 1e-12);
    EXPECT_NEAR(implied[1][1], 0.25, 1e-12);
}
