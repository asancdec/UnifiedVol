// SPDX-License-Identifier: Apache-2.0

#include "Optimization/Helpers.hpp"
#include "Base/Errors/Errors.hpp"
#include "Base/Utils/Detail/Log.hpp"

#include <gtest/gtest.h>
#include <span>
#include <vector>

namespace
{
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
} // namespace

TEST(UnitOptimizationHelpers, ClampBoundsMovesInitialGuessInsideBox)
{
    const SilenceConsoleLog quiet;
    std::vector<double> x{-1.0, 0.5, 3.0};
    const std::vector<double> lower{0.0, 0.0, 0.0};
    const std::vector<double> upper{1.0, 1.0, 2.0};

    uv::opt::clampBounds(
        std::span<double>{x},
        std::span<const double>{lower},
        std::span<const double>{upper}
    );

    EXPECT_DOUBLE_EQ(x[0], 0.0);
    EXPECT_DOUBLE_EQ(x[1], 0.5);
    EXPECT_DOUBLE_EQ(x[2], 2.0);
}

TEST(UnitOptimizationHelpers, ClampLowerAndUpperBoundsIndependently)
{
    const SilenceConsoleLog quiet;
    std::vector<double> x{-1.0, 0.5, 3.0};
    const std::vector<double> lower{0.0, 0.0, 0.0};
    const std::vector<double> upper{0.25, 0.25, 2.0};

    uv::opt::clampLowerBounds(std::span<double>{x}, std::span<const double>{lower});
    uv::opt::clampUpperBounds(std::span<double>{x}, std::span<const double>{upper});

    EXPECT_DOUBLE_EQ(x[0], 0.0);
    EXPECT_DOUBLE_EQ(x[1], 0.25);
    EXPECT_DOUBLE_EQ(x[2], 2.0);
}

TEST(UnitOptimizationHelpers, RejectsInvalidBoundsSpec)
{
    const std::vector<double> lower{1.0};
    const std::vector<double> upper{0.0};

    EXPECT_THROW(
        uv::opt::validateBoundsSpec(
            1,
            std::span<const double>{lower},
            std::span<const double>{upper}
        ),
        uv::errors::UnifiedVolError
    );
}
