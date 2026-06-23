// SPDX-License-Identifier: Apache-2.0

#include "Optimization/Ceres/Optimizer.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include "Base/Errors/Errors.hpp"
#include "Base/Utils/Detail/Log.hpp"

#include <array>
#include <ceres/ceres.h>
#include <gtest/gtest.h>
#include <memory>
#include <span>

namespace uv::tests::integration::optimizer::detail
{
double nloptQuadratic(unsigned n, const double* x, double* grad, void*) noexcept
{
    if (grad)
    {
        grad[0] = 2.0 * (x[0] - 1.5);
        grad[1] = 2.0 * (x[1] + 0.25);
    }

    return (n == 2U) ? (x[0] - 1.5) * (x[0] - 1.5) + (x[1] + 0.25) * (x[1] + 0.25) : 0.0;
}

struct CeresQuadraticResidual
{
    template <typename T> bool operator()(const T* const x, T* residual) const
    {
        residual[0] = x[0] - T{1.5};
        residual[1] = x[1] + T{0.25};
        return true;
    }
};

class SilenceConsoleLog
{
  public:
    SilenceConsoleLog() noexcept
    {
        uv::utils::Log::instance().enableConsole(false);
    }

    ~SilenceConsoleLog() noexcept
    {
        uv::utils::Log::instance().enableConsole(true);
    }
};
} // namespace uv::tests::integration::optimizer::detail

using namespace uv::tests::integration::optimizer::detail;

TEST(IntegrationOptimizerToyProblems, NLoptSolvesBoundedQuadratic)
{
    uv::opt::nlopt::Optimizer<2, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        {.tol = 1e-12,
         .ftolRel = 1e-12,
         .maxEval = 100,
         .verbose = false,
         .paramNames = {"x", "y"}}
    };

    optimizer.setGuessBounds(
        std::array<double, 2>{0.0, 0.0},
        std::array<double, 2>{-10.0, -10.0},
        std::array<double, 2>{10.0, 10.0}
    );
    optimizer.setMinObjective(&nloptQuadratic, nullptr);

    const auto x = optimizer.optimize();

    ASSERT_EQ(x.size(), 2U);
    EXPECT_NEAR(x[0], 1.5, 1e-8);
    EXPECT_NEAR(x[1], -0.25, 1e-8);
}

TEST(IntegrationOptimizerToyProblems, NLoptClampsInitialGuessAndStoresUserValue)
{
    const SilenceConsoleLog quiet;
    uv::opt::nlopt::Optimizer<2, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        {.tol = 1e-12,
         .ftolRel = 1e-12,
         .maxEval = 100,
         .verbose = false,
         .paramNames = {"x", "y"}}
    };

    EXPECT_THROW(optimizer.userValue(), uv::errors::UnifiedVolError);

    optimizer.setUserValue(42.0);
    EXPECT_DOUBLE_EQ(optimizer.userValue(), 42.0);

    optimizer.setGuessBounds(
        std::array<double, 2>{100.0, -100.0},
        std::array<double, 2>{-1.0, -1.0},
        std::array<double, 2>{1.0, 1.0}
    );
    optimizer.setMinObjective(&nloptQuadratic, nullptr);

    const auto x = optimizer.optimize();

    ASSERT_EQ(x.size(), 2U);
    EXPECT_NEAR(x[0], 1.0, 1e-10);
    EXPECT_NEAR(x[1], -0.25, 1e-8);
}

TEST(IntegrationOptimizerToyProblems, NLoptRejectsInvalidBounds)
{
    uv::opt::nlopt::Optimizer<2, uv::opt::nlopt::Algorithm::LD_SLSQP> optimizer{
        {.tol = 1e-12,
         .ftolRel = 1e-12,
         .maxEval = 100,
         .verbose = false,
         .paramNames = {"x", "y"}}
    };

    EXPECT_THROW(
        optimizer.setGuessBounds(
            std::array<double, 2>{0.0, 0.0},
            std::array<double, 2>{-1.0, 1.0},
            std::array<double, 2>{1.0, -1.0}
        ),
        uv::errors::UnifiedVolError
    );
}

TEST(IntegrationOptimizerToyProblems, CeresSolvesBoundedQuadratic)
{
    uv::opt::ceres::Optimizer optimizer{
        {.maxEval = 50,
         .functionTol = 1e-12,
         .paramTol = 1e-12,
         .gradientTol = 1e-12,
         .lossScale = 1.0,
         .paramNames = {"x", "y"},
         .verbosity = uv::opt::ceres::Verbosity::None,
         .numThreads = 1}
    };
    const std::array<double, 2> init{0.0, 0.0};
    const std::array<double, 2> lower{-10.0, -10.0};
    const std::array<double, 2> upper{10.0, 10.0};

    optimizer.initialize(
        std::span<const double>{init},
        std::span<const double>{lower},
        std::span<const double>{upper}
    );
    optimizer.beginRun();
    optimizer.addResidualBlock(
        std::make_unique<::ceres::AutoDiffCostFunction<CeresQuadraticResidual, 2, 2>>(
            new CeresQuadraticResidual{}
        )
    );

    const auto x = optimizer.solve();

    ASSERT_EQ(x.size(), 2U);
    EXPECT_NEAR(x[0], 1.5, 1e-10);
    EXPECT_NEAR(x[1], -0.25, 1e-10);
}

TEST(IntegrationOptimizerToyProblems, CeresEnforcesLifecycle)
{
    uv::opt::ceres::Optimizer optimizer{
        {.maxEval = 50,
         .functionTol = 1e-12,
         .paramTol = 1e-12,
         .gradientTol = 1e-12,
         .lossScale = 1.0,
         .paramNames = {"x", "y"},
         .verbosity = uv::opt::ceres::Verbosity::None,
         .numThreads = 1}
    };
    const std::array<double, 2> init{0.0, 0.0};
    const std::array<double, 2> lower{-10.0, -10.0};
    const std::array<double, 2> upper{10.0, 10.0};

    EXPECT_THROW(optimizer.params(), uv::errors::UnifiedVolError);
    EXPECT_THROW(optimizer.beginRun(), uv::errors::UnifiedVolError);

    optimizer.initialize(
        std::span<const double>{init},
        std::span<const double>{lower},
        std::span<const double>{upper}
    );

    EXPECT_THROW(optimizer.params(), uv::errors::UnifiedVolError);
    EXPECT_THROW(
        optimizer.addResidualBlock(
            std::make_unique<::ceres::AutoDiffCostFunction<CeresQuadraticResidual, 2, 2>>(
                new CeresQuadraticResidual{}
            )
        ),
        uv::errors::UnifiedVolError
    );
}
