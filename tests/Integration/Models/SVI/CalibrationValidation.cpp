// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"
#include "Core/Matrix.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <gtest/gtest.h>
#include <limits>
#include <span>
#include <vector>

namespace
{
using SVIOptimizer = uv::opt::nlopt::Optimizer<4, uv::opt::nlopt::Algorithm::LD_SLSQP>;

SVIOptimizer makeOptimizer()
{
    const uv::models::svi::Config config{
        .objectiveTol = 1e-12,
        .maxEval = 100,
        .verbose = false,
        .printParams = false
    };
    return SVIOptimizer{uv::models::svi::detail::makeNLoptConfig(config)};
}
} // namespace

TEST(IntegrationSVICalibrationValidation, RejectsNonIncreasingMaturities)
{
    const std::vector<double> maturities{0.5, 0.5};
    const uv::core::Matrix<double> logKF{2, 3, 0.0};
    const uv::core::Matrix<double> totalVariance{2, 3, 0.04};
    const auto optimizer = makeOptimizer();

    EXPECT_THROW(
        uv::models::svi::calibrate(
            std::span<const double>{maturities},
            logKF,
            totalVariance,
            optimizer
        ),
        uv::errors::UnifiedVolError
    );
}

TEST(IntegrationSVICalibrationValidation, RejectsShapeMismatch)
{
    const std::vector<double> maturities{0.5, 1.0};
    const uv::core::Matrix<double> logKF{2, 3, 0.0};
    const uv::core::Matrix<double> totalVariance{2, 2, 0.04};
    const auto optimizer = makeOptimizer();

    EXPECT_THROW(
        uv::models::svi::calibrate(
            std::span<const double>{maturities},
            logKF,
            totalVariance,
            optimizer
        ),
        uv::errors::UnifiedVolError
    );
}

TEST(IntegrationSVICalibrationValidation, RejectsInvalidSliceData)
{
    const std::vector<double> maturities{0.5};
    uv::core::Matrix<double> logKF{1, 3, 0.0};
    uv::core::Matrix<double> totalVariance{1, 3, 0.04};
    const auto optimizer = makeOptimizer();

    logKF[0][0] = -0.1;
    logKF[0][1] = -0.1;
    logKF[0][2] = 0.1;
    EXPECT_THROW(
        uv::models::svi::calibrate(
            std::span<const double>{maturities},
            logKF,
            totalVariance,
            optimizer
        ),
        uv::errors::UnifiedVolError
    );

    logKF[0][0] = -0.1;
    logKF[0][1] = 0.0;
    logKF[0][2] = 0.1;
    totalVariance[0][1] = -0.01;
    EXPECT_THROW(
        uv::models::svi::calibrate(
            std::span<const double>{maturities},
            logKF,
            totalVariance,
            optimizer
        ),
        uv::errors::UnifiedVolError
    );

    totalVariance[0][1] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(
        uv::models::svi::calibrate(
            std::span<const double>{maturities},
            logKF,
            totalVariance,
            optimizer
        ),
        uv::errors::UnifiedVolError
    );
}
