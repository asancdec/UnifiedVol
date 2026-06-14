// SPDX-License-Identifier: Apache-2.0

#include "../../Golden.hpp"
#include "Core/Matrix.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <span>
#include <vector>

namespace
{
struct ErrorDiagnostics
{
    double meanAbs{};
    double maxAbs{};
    double rmse{};
};

struct SyntheticSVICase
{
    std::vector<double> maturities{0.5, 1.0};
    std::vector<double> logKF{-0.25, -0.10, 0.0, 0.10, 0.25};
    uv::core::Matrix<double> logKFMatrix{maturities.size(), logKF.size()};
    uv::core::Matrix<double> totalVariance{maturities.size(), logKF.size()};
    uv::Vector<uv::models::svi::Params<double>> truth{
        {0.5, 0.025, 0.30, -0.35, 0.02, 0.30},
        {1.0, 0.040, 0.35, -0.25, 0.04, 0.35}
    };
};

SyntheticSVICase makeSyntheticSVICase()
{
    SyntheticSVICase data;

    for (std::size_t i = 0; i < data.maturities.size(); ++i)
    {
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            data.logKFMatrix[i][j] = data.logKF[j];
            const auto& p = data.truth[i];
            data.totalVariance[i][j] = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
        }
    }

    return data;
}

uv::Vector<uv::models::svi::Params<double>>
calibrateSynthetic(const SyntheticSVICase& data)
{
    const uv::models::svi::Config config{
        .objectiveTol = 1e-12,
        .maxEval = 500,
        .verbose = false,
        .printParams = false
    };
    const uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        uv::models::svi::detail::makeNLoptConfig(config)
    };

    return uv::models::svi::calibrate(
        std::span<const double>{data.maturities},
        data.logKFMatrix,
        data.totalVariance,
        optimizer
    );
}

void expectParamsNear(
    const uv::models::svi::Params<double>& actual,
    const uv::models::svi::Params<double>& expected,
    const double tolerance
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

ErrorDiagnostics totalVarianceDiagnostics(
    const uv::Vector<uv::models::svi::Params<double>>& calibrated,
    const SyntheticSVICase& data
)
{
    double sumAbs{};
    double sumSquared{};
    double maxAbs{};
    std::size_t count{};

    for (std::size_t i = 0; i < calibrated.size(); ++i)
    {
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            const auto& p = calibrated[i];
            const double fitted = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
            const double error = std::abs(fitted - data.totalVariance[i][j]);

            sumAbs += error;
            sumSquared += error * error;
            maxAbs = std::max(maxAbs, error);
            ++count;
        }
    }

    return {
        sumAbs / static_cast<double>(count),
        maxAbs,
        std::sqrt(sumSquared / static_cast<double>(count))
    };
}
} // namespace

TEST(RegressionSVICalibration, RecoversSyntheticSurfaceShapeAndVols)
{
    const auto data = makeSyntheticSVICase();
    const auto calibrated = calibrateSynthetic(data);
    const auto golden = uv::tests::golden::readSyntheticSVICalibration(
        "tests/Golden/synthetic_svi_calibration.json"
    );

    ASSERT_EQ(calibrated.size(), data.truth.size());
    ASSERT_EQ(calibrated.size(), golden.calibratedParams.size());
    for (std::size_t i = 0; i < calibrated.size(); ++i)
        expectParamsNear(calibrated[i], golden.calibratedParams[i], golden.tolerance);

    const ErrorDiagnostics errors{totalVarianceDiagnostics(calibrated, data)};
    EXPECT_NEAR(errors.meanAbs, golden.meanAbsTotalVarianceError, golden.tolerance)
        << "metric=meanAbsTotalVariance";
    EXPECT_NEAR(errors.maxAbs, golden.maxAbsTotalVarianceError, golden.tolerance)
        << "metric=maxAbsTotalVariance";
    EXPECT_NEAR(errors.rmse, golden.rmseTotalVarianceError, golden.tolerance)
        << "metric=rmseTotalVariance";

    for (std::size_t i = 0; i < calibrated.size(); ++i)
    {
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            const auto& p = calibrated[i];
            const double fitted = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
            EXPECT_NEAR(fitted, data.totalVariance[i][j], golden.tolerance)
                << "slice=" << i << " strike=" << j << " logKF=" << data.logKF[j]
                << " actual=" << fitted << " expected=" << data.totalVariance[i][j];
        }
    }
}

TEST(RegressionSVIStaticArbitrage, CalibratedParamsSatisfyAllConstraints)
{
    constexpr double eps = 1e-12;
    constexpr double wingDelta = 0.15;

    const auto data = makeSyntheticSVICase();
    const auto calibrated = calibrateSynthetic(data);
    const auto [minIt, maxIt] = std::minmax_element(data.logKF.begin(), data.logKF.end());

    std::vector<double> constraintLogKF = data.logKF;
    constraintLogKF.push_back(*minIt - wingDelta);
    constraintLogKF.push_back(*maxIt + wingDelta);

    ASSERT_EQ(calibrated.size(), data.truth.size());

    for (std::size_t i = 0; i < calibrated.size(); ++i)
    {
        const auto& p = calibrated[i];
        const double sqrtOneMinusRho2 = std::sqrt(1.0 - p.rho * p.rho);
        const double wMin = p.a + p.b * p.sigma * sqrtOneMinusRho2;
        const double leftSlope = p.b * (1.0 + p.rho);
        const double rightSlope = p.b * (1.0 - p.rho);

        EXPECT_GT(wMin, eps) << "slice " << i;
        EXPECT_GT(leftSlope, eps) << "slice " << i;
        EXPECT_GT(rightSlope, eps) << "slice " << i;
        EXPECT_LT(leftSlope, 2.0) << "slice " << i;
        EXPECT_LT(rightSlope, 2.0) << "slice " << i;

        for (double k : constraintLogKF)
        {
            const double w =
                uv::models::svi::totalVariance(p.a, p.b, p.rho, p.m, p.sigma, k);
            const double g = uv::models::svi::gk(p.a, p.b, p.rho, p.m, p.sigma, k);

            EXPECT_TRUE(std::isfinite(w)) << "slice " << i << " k=" << k;
            EXPECT_TRUE(std::isfinite(g)) << "slice " << i << " k=" << k;
            EXPECT_GT(w, eps) << "slice " << i << " k=" << k;
            EXPECT_GE(g, -1e-10) << "slice " << i << " k=" << k;

            if (i > 0)
            {
                const auto& prev = calibrated[i - 1];
                const double prevW = uv::models::svi::totalVariance(
                    prev.a,
                    prev.b,
                    prev.rho,
                    prev.m,
                    prev.sigma,
                    k
                );

                EXPECT_GE(w + 1e-10, prevW) << "slice " << i << " k=" << k;
            }
        }
    }
}
