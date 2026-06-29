// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/Config.hpp"
#include "Models/SVI/Math.hpp"
#include "Optimization/NLopt/Optimizer.hpp"
#include "Support/Assertions.hpp"
#include "Support/Golden.hpp"
#include "Support/Models/SVI.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <span>
#include <vector>

namespace uv::tests::regression::svi_synthetic::detail
{
namespace assertions = uv::tests::assertions;
namespace svi_support = uv::tests::models::svi;

uv::Vector<uv::models::svi::Params<double>>
calibrateSynthetic(const svi_support::SyntheticSurfaceCase& data)
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

} // namespace uv::tests::regression::svi_synthetic::detail

using namespace uv::tests::regression::svi_synthetic::detail;

TEST(RegressionSVICalibration, RecoversSyntheticSurfaceShapeAndVols)
{
    const auto data = svi_support::makeDefaultSyntheticSurfaceCase();
    const auto calibrated = calibrateSynthetic(data);
    const auto golden = uv::tests::golden::readSyntheticSVICalibration(
        "tests/Golden/synthetic_svi_calibration.json"
    );

    ASSERT_EQ(calibrated.size(), data.truth.size());
    ASSERT_EQ(calibrated.size(), golden.calibratedParams.size());
    for (std::size_t i = 0; i < calibrated.size(); ++i)
        assertions::expectSVIParamsNear(
            calibrated[i],
            golden.calibratedParams[i],
            golden.tolerance
        );

    const uv::tests::ErrorDiagnostics errors{
        svi_support::totalVarianceDiagnostics(calibrated, data)
    };
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

    const auto data = svi_support::makeDefaultSyntheticSurfaceCase();
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
