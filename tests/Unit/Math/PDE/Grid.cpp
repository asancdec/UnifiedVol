// SPDX-License-Identifier: Apache-2.0

#include "Math/PDE/Grid.hpp"
#include "Base/Errors/Errors.hpp"

#include <array>
#include <gtest/gtest.h>

TEST(UnitMathPDEGrid, StoresCoordinatesAndSteps)
{
    const std::array<double, 4> x{-1.0, -0.25, 0.5, 2.0};
    const uv::math::pde::Grid<double, 4> grid{x};

    EXPECT_DOUBLE_EQ(grid.x()[0], -1.0);
    EXPECT_DOUBLE_EQ(grid.x()[3], 2.0);
    EXPECT_DOUBLE_EQ(grid.dx()[0], 0.75);
    EXPECT_DOUBLE_EQ(grid.dx()[2], 1.5);
}

TEST(UnitMathPDEGrid, GeneratesUniformGrid)
{
    const auto grid = uv::math::pde::generateUniformGrid<double, 5>(-1.0, 1.0);

    EXPECT_DOUBLE_EQ(grid.x()[0], -1.0);
    EXPECT_DOUBLE_EQ(grid.x()[2], 0.0);
    EXPECT_DOUBLE_EQ(grid.x()[4], 1.0);
    EXPECT_DOUBLE_EQ(grid.dx()[0], 0.5);
    EXPECT_DOUBLE_EQ(grid.dx()[3], 0.5);
}

TEST(UnitMathPDEGrid, GeneratesCenteredSinhGrid)
{
    const auto grid = uv::math::pde::generateCenteredSinHGrid<double, 5>(-2.0, 2.0, 1.0);

    EXPECT_NEAR(grid.x()[0], -2.0, 1e-15);
    EXPECT_NEAR(grid.x()[2], 0.0, 1e-15);
    EXPECT_NEAR(grid.x()[4], 2.0, 1e-15);
    EXPECT_LT(grid.dx()[1], grid.dx()[0]);
    EXPECT_NEAR(grid.dx()[1], grid.dx()[2], 1e-15);
    EXPECT_LT(grid.dx()[2], grid.dx()[3]);
}

TEST(UnitMathPDEGrid, RejectsNonMonotonicInput)
{
    const std::array<double, 4> x{0.0, 1.0, 0.5, 2.0};

    EXPECT_THROW((uv::math::pde::Grid<double, 4>{x}), uv::errors::UnifiedVolError);
}

TEST(UnitMathPDEGrid, CenteredSinhGridWithTinyBetaFallsBackToUniformGrid)
{
    const auto grid = uv::math::pde::generateCenteredSinHGrid<double, 5>(-1.0, 1.0, 0.0);

    EXPECT_DOUBLE_EQ(grid.x()[0], -1.0);
    EXPECT_DOUBLE_EQ(grid.x()[2], 0.0);
    EXPECT_DOUBLE_EQ(grid.x()[4], 1.0);
    EXPECT_DOUBLE_EQ(grid.dx()[0], 0.5);
    EXPECT_DOUBLE_EQ(grid.dx()[3], 0.5);
}

TEST(UnitMathPDEGrid, RejectsInvalidCenteredSinhInputs)
{
    EXPECT_THROW(
        (uv::math::pde::generateCenteredSinHGrid<double, 5>(-1.0, 2.0, 1.0)),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        (uv::math::pde::generateCenteredSinHGrid<double, 5>(-1.0, 1.0, -0.1)),
        uv::errors::UnifiedVolError
    );
}
