// SPDX-License-Identifier: Apache-2.0

#include "Optimization/Cost.hpp"
#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(UnitOptimizationCost, WeightsATMEmphasizesAtTheMoneyPoint)
{
    const std::vector<double> logKF{-0.2, 0.0, 0.2};
    std::vector<double> weights(3);

    uv::opt::cost::weightsATM<double>(
        std::span<const double>{logKF},
        {.wATM = 4.0, .k0 = 0.1},
        std::span<double>{weights}
    );

    EXPECT_LT(weights[0], weights[1]);
    EXPECT_DOUBLE_EQ(weights[1], 2.0);
    EXPECT_LT(weights[2], weights[1]);
}

TEST(UnitOptimizationCost, RejectsWeightBelowOne)
{
    const std::vector<double> logKF{0.0};
    std::vector<double> weights(1);

    EXPECT_THROW(
        uv::opt::cost::weightsATM<double>(
            std::span<const double>{logKF},
            {.wATM = 0.5, .k0 = 0.1},
            std::span<double>{weights}
        ),
        uv::errors::UnifiedVolError
    );
}
