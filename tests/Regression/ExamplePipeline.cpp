// SPDX-License-Identifier: Apache-2.0

#include "Golden.hpp"
#include "IO/Load.hpp"
#include "Models/Heston/BuildSurface.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Optimization/NLopt/Optimizer.hpp"
#include "Support/Tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <gtest/gtest.h>

namespace
{
struct ErrorDiagnostics
{
    double meanAbs{};
    double maxAbs{};
    double rmse{};
};

template <std::floating_point T>
void expectFiniteNonNegativeVolSurface(const uv::core::VolSurface<T>& surface)
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

template <std::floating_point T> ErrorDiagnostics volErrorDiagnostics(
    const uv::core::VolSurface<T>& lhs,
    const uv::core::VolSurface<T>& rhs
)
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
    const uv::models::svi::Params<T>& actual,
    const uv::models::svi::Params<T>& expected,
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
    const uv::models::heston::Params<T>& actual,
    const uv::models::heston::Params<T>& expected,
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
    const uv::core::VolSurface<T>& surface,
    const uv::Vector<uv::tests::golden::VolPoint>& expected,
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
    const uv::core::Matrix<T>& prices,
    const uv::core::VolSurface<T>& surface,
    const uv::core::Curve<T>& curve
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
            EXPECT_GE(price + uv::tests::tolerance::NoArb, intrinsic)
                << "row=" << i << " col=" << j;
            EXPECT_LE(
                price,
                discountFactors[i] * forwards[i] + uv::tests::tolerance::NoArb
            ) << "row="
              << i << " col=" << j;

            if (j > 0)
            {
                EXPECT_LE(price, prices[i][j - 1] + uv::tests::tolerance::NoArb)
                    << "row=" << i << " col=" << j;
            }
            if (j > 0 && j + 1 < prices.cols())
            {
                const T leftSlope =
                    (prices[i][j] - prices[i][j - 1]) / (strikes[j] - strikes[j - 1]);
                const T rightSlope =
                    (prices[i][j + 1] - prices[i][j]) / (strikes[j + 1] - strikes[j]);
                EXPECT_GE(rightSlope + uv::tests::tolerance::NoArb, leftSlope)
                    << "row=" << i << " col=" << j;
            }
        }
    }
}
} // namespace

TEST(RegressionExamplePipeline, MainCppResultsRemainStable)
{
    using Real = double;

    const std::filesystem::path path{"data/VolSurface_SPY_04072011.csv"};
    const auto golden =
        uv::tests::golden::readExamplePipeline("tests/Golden/example_pipeline.json");
    const Real tolerance{golden.tolerance};
    const uv::core::MarketData<Real> marketData{
        .interestRate = 0.0,
        .dividendYield = 0.0,
        .spot = 485.77548
    };
    const uv::models::svi::Config sviConfig{
        .objectiveTol = 1e-12,
        .maxEval = 10000,
        .verbose = false,
        .printParams = false
    };
    const uv::models::heston::calibrate::Config hestonConfig{
        .tolerance = 1e-11,
        .maxEval = 10000,
        .verbosity = uv::opt::ceres::Verbosity::None
    };

    const auto marketState = uv::io::load::marketState(path, marketData);

    const uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP> sviOptimizer{
        uv::models::svi::detail::makeNLoptConfig(sviConfig)
    };
    const auto sviParams =
        uv::models::svi::calibrate(marketState.volSurface, sviOptimizer, false);
    const auto sviVolSurface =
        uv::models::svi::buildSurface(marketState.volSurface, sviParams);

    const auto hestonParams = uv::models::heston::calibrate::calibrate<Real>(
        sviVolSurface,
        marketState.interestCurve,
        hestonConfig
    );

    uv::models::heston::price::Pricer<Real> hestonPricer{};
    hestonPricer.setParams(hestonParams);
    const auto hestonVolSurface = uv::models::heston::buildSurface(
        sviVolSurface,
        marketState.interestCurve,
        hestonPricer
    );

    ASSERT_EQ(marketState.volSurface.numMaturities(), 11U);
    ASSERT_EQ(marketState.volSurface.numStrikes(), 17U);
    ASSERT_EQ(sviVolSurface.numMaturities(), marketState.volSurface.numMaturities());
    ASSERT_EQ(sviVolSurface.numStrikes(), marketState.volSurface.numStrikes());
    ASSERT_EQ(hestonVolSurface.numMaturities(), marketState.volSurface.numMaturities());
    ASSERT_EQ(hestonVolSurface.numStrikes(), marketState.volSurface.numStrikes());

    expectFiniteNonNegativeVolSurface(marketState.volSurface);
    expectFiniteNonNegativeVolSurface(sviVolSurface);
    expectFiniteNonNegativeVolSurface(hestonVolSurface);

    ASSERT_EQ(sviParams.size(), golden.sviParams.size());
    for (std::size_t i = 0; i < sviParams.size(); ++i)
        expectSVIParamsNear(sviParams[i], golden.sviParams[i], tolerance);
    expectHestonParamsNear(hestonParams, golden.hestonParams, tolerance);
    expectVolPointsNear(marketState.volSurface, golden.marketVols, tolerance);
    expectVolPointsNear(sviVolSurface, golden.sviVols, tolerance);
    expectVolPointsNear(hestonVolSurface, golden.hestonVols, tolerance);
    expectCallSurfaceNoArbitrage(
        hestonPricer.callPrice(sviVolSurface, marketState.interestCurve),
        sviVolSurface,
        marketState.interestCurve
    );

    const ErrorDiagnostics errors{volErrorDiagnostics(hestonVolSurface, sviVolSurface)};

    EXPECT_NEAR(errors.meanAbs, golden.meanAbsVolError, tolerance);
    EXPECT_NEAR(errors.maxAbs, golden.maxAbsVolError, tolerance);
    EXPECT_NEAR(errors.rmse, golden.rmseVolError, tolerance);
}
