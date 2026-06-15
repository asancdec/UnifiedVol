// SPDX-License-Identifier: Apache-2.0

#include "Optimization/Helpers.hpp"
#include "Base/Errors/Errors.hpp"
#include "Base/Utils/Detail/Log.hpp"

#include <gtest/gtest.h>
#include <limits>
#include <span>
#include <string_view>
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

TEST(UnitOptimizationHelpers, ValidatesFullBoundsInputs)
{
    const std::vector<double> x{0.25, 0.5};
    const std::vector<double> lower{0.0, 0.0};
    const std::vector<double> upper{1.0, 1.0};

    EXPECT_NO_THROW(uv::opt::validateBounds(x, lower, upper));
    EXPECT_NO_THROW(uv::opt::validateBoundsSpec(2, lower, upper));
}

TEST(UnitOptimizationHelpers, RejectsInvalidFullBoundsInputs)
{
    const std::vector<double> x{0.25, 0.5};
    const std::vector<double> lower{0.0};
    const std::vector<double> upper{1.0, 1.0};
    const std::vector<double> nonFinite{0.25, std::numeric_limits<double>::quiet_NaN()};

    EXPECT_THROW(uv::opt::validateBounds({}, {}, {}), uv::errors::UnifiedVolError);
    EXPECT_THROW(
        uv::opt::validateBounds(nonFinite, upper, upper),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(uv::opt::validateBounds(x, lower, upper), uv::errors::UnifiedVolError);
}

TEST(UnitOptimizationHelpers, ValidatesOneSidedBounds)
{
    const std::vector<double> x{0.25, 0.5};
    const std::vector<double> lower{0.0, 0.0};
    const std::vector<double> upper{1.0, 1.0};

    EXPECT_NO_THROW(uv::opt::validateLowerBounds(x, lower));
    EXPECT_NO_THROW(uv::opt::validateLowerBoundsSpec(2, lower));
    EXPECT_NO_THROW(uv::opt::validateUpperBounds(x, upper));
    EXPECT_NO_THROW(uv::opt::validateUpperBoundsSpec(2, upper));
}

TEST(UnitOptimizationHelpers, RejectsInvalidOneSidedBounds)
{
    const std::vector<double> x{0.25, 0.5};
    const std::vector<double> shortBounds{0.0};
    const std::vector<double> nonFinite{0.25, std::numeric_limits<double>::infinity()};

    EXPECT_THROW(
        uv::opt::validateLowerBounds({}, shortBounds),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::opt::validateLowerBounds(nonFinite, nonFinite),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::opt::validateLowerBounds(x, shortBounds),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::opt::validateUpperBounds(x, shortBounds),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitOptimizationHelpers, WarnBoundsHitAcceptsOptionalBounds)
{
    const SilenceConsoleLog quiet;
    const std::vector<double> x{0.0, 1.0};
    const std::vector<double> lower{0.0, -1.0};
    const std::vector<double> upper{1.0, 1.0};

    EXPECT_NO_THROW(uv::opt::warnBoundsHit(x, std::nullopt, std::nullopt));
    EXPECT_NO_THROW(uv::opt::warnBoundsHit(x, lower, std::nullopt));
    EXPECT_NO_THROW(uv::opt::warnBoundsHit(x, std::nullopt, upper));
    EXPECT_NO_THROW(uv::opt::warnBoundsHit(
        std::span<const double>{x},
        std::span<const double>{lower},
        std::span<const double>{upper}
    ));
}

TEST(UnitOptimizationHelpers, LogResultsHandlesNamedAndUnnamedParameters)
{
    const SilenceConsoleLog quiet;
    const std::vector<double> x{1.0, 2.0};
    const std::vector<std::string_view> names{"a", "b"};

    EXPECT_NO_THROW(uv::opt::logResults(x, {}, 0.25, 3U, 1.5, true));
    EXPECT_NO_THROW(uv::opt::logResults(x, names, 0.25, 3U, 1.5, false, "FAILURE"));
}

TEST(UnitOptimizationHelpers, LogResultsRejectsMismatchedParamNames)
{
    const SilenceConsoleLog quiet;
    const std::vector<double> x{1.0, 2.0};
    const std::vector<std::string_view> names{"a"};

    EXPECT_THROW(
        uv::opt::logResults(x, names, 0.25, 3U, 1.5, true),
        uv::errors::UnifiedVolError
    );
}
