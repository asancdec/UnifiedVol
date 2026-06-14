// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/SVI/Params.hpp"
#include "Support/Diagnostics.hpp"
#include "Support/Tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <gtest/gtest.h>

namespace uv::tests::assertions
{
template <std::floating_point T>
void expectFiniteNonNegativeVolSurface(const core::VolSurface<T>& surface)
{
    for (std::size_t i = 0; i < surface.numMaturities(); ++i)
    {
        for (std::size_t j = 0; j < surface.numStrikes(); ++j)
        {
            EXPECT_TRUE(std::isfinite(surface.vol()[i][j]));
            EXPECT_GE(surface.vol()[i][j], T{0});
        }
    }
}

template <std::floating_point T> ErrorDiagnostics
volErrorDiagnostics(const core::VolSurface<T>& lhs, const core::VolSurface<T>& rhs)
{
    T sum{0};
    T sumSquared{0};
    T maxError{0};
    std::size_t count{0};

    for (std::size_t i = 0; i < lhs.numMaturities(); ++i)
    {
        for (std::size_t j = 0; j < lhs.numStrikes(); ++j)
        {
            const T error = std::abs(lhs.vol()[i][j] - rhs.vol()[i][j]);
            sum += error;
            sumSquared += error * error;
            maxError = std::max(maxError, error);
            ++count;
        }
    }

    return {
        static_cast<double>(sum / static_cast<T>(count)),
        static_cast<double>(maxError),
        static_cast<double>(std::sqrt(sumSquared / static_cast<T>(count)))
    };
}

template <std::floating_point T> void expectSVIParamsNear(
    const models::svi::Params<T>& actual,
    const models::svi::Params<T>& expected,
    const T tolerance
)
{
    EXPECT_NEAR(actual.t, expected.t, tolerance) << "field=t";
    EXPECT_NEAR(actual.a, expected.a, tolerance) << "field=a t=" << expected.t;
    EXPECT_NEAR(actual.b, expected.b, tolerance) << "field=b t=" << expected.t;
    EXPECT_NEAR(actual.rho, expected.rho, tolerance) << "field=rho t=" << expected.t;
    EXPECT_NEAR(actual.m, expected.m, tolerance) << "field=m t=" << expected.t;
    EXPECT_NEAR(actual.sigma, expected.sigma, tolerance)
        << "field=sigma t=" << expected.t;
}

template <std::floating_point T> void expectHestonParamsNear(
    const models::heston::Params<T>& actual,
    const models::heston::Params<T>& expected,
    const T tolerance
)
{
    EXPECT_NEAR(actual.kappa, expected.kappa, tolerance) << "field=kappa";
    EXPECT_NEAR(actual.theta, expected.theta, tolerance) << "field=theta";
    EXPECT_NEAR(actual.sigma, expected.sigma, tolerance) << "field=sigma";
    EXPECT_NEAR(actual.rho, expected.rho, tolerance) << "field=rho";
    EXPECT_NEAR(actual.v0, expected.v0, tolerance) << "field=v0";
}

template <std::floating_point T> void expectVolPointsNear(
    const core::VolSurface<T>& surface,
    const auto& expected,
    const T tolerance
)
{
    for (const auto& point : expected)
    {
        ASSERT_LT(point.maturity, surface.numMaturities());
        ASSERT_LT(point.strike, surface.numStrikes());
        EXPECT_NEAR(surface.vol()[point.maturity][point.strike], point.value, tolerance)
            << "maturity=" << point.maturity << " strike=" << point.strike;
    }
}

template <std::floating_point T> void expectCallSurfaceNoArbitrage(
    const core::Matrix<T>& prices,
    const core::VolSurface<T>& surface,
    const core::Curve<T>& curve
)
{
    const auto discountFactors = curve.interpolateDF(surface.maturities());
    const auto forwards = surface.forwards();
    const auto strikes = surface.strikes();

    ASSERT_EQ(prices.rows(), surface.numMaturities());
    ASSERT_EQ(prices.cols(), surface.numStrikes());

    for (std::size_t i = 0; i < prices.rows(); ++i)
    {
        for (std::size_t j = 0; j < prices.cols(); ++j)
        {
            const T price = prices[i][j];
            const T intrinsic =
                discountFactors[i] * std::max(forwards[i] - strikes[j], T{0});
            EXPECT_TRUE(std::isfinite(price)) << "row=" << i << " col=" << j;
            EXPECT_GE(price + tolerance::NoArb, intrinsic) << "row=" << i << " col=" << j;
            EXPECT_LE(price, discountFactors[i] * forwards[i] + tolerance::NoArb)
                << "row=" << i << " col=" << j;

            if (j > 0)
            {
                EXPECT_LE(price, prices[i][j - 1] + tolerance::NoArb)
                    << "row=" << i << " col=" << j;
            }
            if (j > 0 && j + 1 < prices.cols())
            {
                const T leftSlope =
                    (prices[i][j] - prices[i][j - 1]) / (strikes[j] - strikes[j - 1]);
                const T rightSlope =
                    (prices[i][j + 1] - prices[i][j]) / (strikes[j + 1] - strikes[j]);
                EXPECT_GE(rightSlope + tolerance::NoArb, leftSlope)
                    << "row=" << i << " col=" << j;
            }
        }
    }
}
} // namespace uv::tests::assertions
