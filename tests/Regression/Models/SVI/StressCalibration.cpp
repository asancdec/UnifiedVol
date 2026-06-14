// SPDX-License-Identifier: Apache-2.0

#include "../../../Support/Tolerances.hpp"
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
struct StressSVICase
{
    std::vector<double> maturities;
    std::vector<double> logKF;
    uv::Vector<uv::models::svi::Params<double>> truth;
};

uv::core::Matrix<double> makeLogKFMatrix(const StressSVICase& data)
{
    uv::core::Matrix<double> logKF{data.maturities.size(), data.logKF.size()};
    for (std::size_t i = 0; i < logKF.rows(); ++i)
    {
        for (std::size_t j = 0; j < logKF.cols(); ++j)
            logKF[i][j] = data.logKF[j];
    }
    return logKF;
}

uv::core::Matrix<double> makeTotalVariance(const StressSVICase& data)
{
    uv::core::Matrix<double> totalVariance{data.maturities.size(), data.logKF.size()};
    for (std::size_t i = 0; i < totalVariance.rows(); ++i)
    {
        for (std::size_t j = 0; j < totalVariance.cols(); ++j)
        {
            const auto& p = data.truth[i];
            totalVariance[i][j] = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
        }
    }
    return totalVariance;
}

uv::Vector<uv::models::svi::Params<double>> calibrateStressCase(
    const StressSVICase& data,
    const uv::core::Matrix<double>& logKF,
    const uv::core::Matrix<double>& totalVariance
)
{
    const uv::models::svi::Config config{
        .objectiveTol = 1e-12,
        .maxEval = 2000,
        .verbose = false,
        .printParams = false
    };
    const uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        uv::models::svi::detail::makeNLoptConfig(config)
    };

    return uv::models::svi::calibrate(
        std::span<const double>{data.maturities},
        logKF,
        totalVariance,
        optimizer
    );
}

void expectStressFitWithinTolerance(const StressSVICase& data, const double tolerance)
{
    const auto logKF = makeLogKFMatrix(data);
    const auto totalVariance = makeTotalVariance(data);
    const auto calibrated = calibrateStressCase(data, logKF, totalVariance);

    ASSERT_EQ(calibrated.size(), data.truth.size());
    for (std::size_t i = 0; i < calibrated.size(); ++i)
    {
        const auto& p = calibrated[i];
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            const double fitted = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
            EXPECT_NEAR(fitted, totalVariance[i][j], tolerance)
                << "slice=" << i << " strike=" << j;
        }
    }
}
} // namespace

TEST(StressSVICalibration, RecoversLowVolAndSteepSkewSyntheticSurfaces)
{
    const std::vector<StressSVICase> cases{
        {{0.05, 0.5, 2.0},
         {-0.35, -0.15, 0.0, 0.15, 0.35},
         {{0.05, 0.0005, 0.010, -0.10, 0.00, 0.050},
          {0.5, 0.0040, 0.020, -0.10, 0.00, 0.100},
          {2.0, 0.0200, 0.030, -0.05, 0.00, 0.160}}},
        {{0.08, 0.25, 1.0},
         {-0.60, -0.30, -0.10, 0.0, 0.15, 0.35, 0.60},
         {{0.08, 0.010, 0.22, -0.80, -0.05, 0.22},
          {0.25, 0.030, 0.30, -0.75, -0.04, 0.28},
          {1.0, 0.080, 0.36, -0.65, -0.02, 0.36}}}
    };

    for (const auto& data : cases)
        expectStressFitWithinTolerance(data, 5e-5);
}
