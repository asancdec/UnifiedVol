// SPDX-License-Identifier: Apache-2.0

#include "Core/Matrix.hpp"

#include <gtest/gtest.h>

TEST(CoreMatrix, ConstructsWithShapeAndFillValue)
{
    uv::core::Matrix<double> m{2, 3, 1.5};

    EXPECT_FALSE(m.empty());
    EXPECT_EQ(m.rows(), 2U);
    EXPECT_EQ(m.cols(), 3U);

    for (std::size_t i = 0; i < m.rows(); ++i)
    {
        for (std::size_t j = 0; j < m.cols(); ++j)
        {
            EXPECT_DOUBLE_EQ(m[i][j], 1.5);
        }
    }
}

TEST(CoreMatrix, SupportsElementwiseMatrixArithmetic)
{
    uv::core::Matrix<double> lhs{2, 2};
    lhs[0][0] = 1.0;
    lhs[0][1] = 2.0;
    lhs[1][0] = 3.0;
    lhs[1][1] = 4.0;

    uv::core::Matrix<double> rhs{2, 2, 0.5};

    const auto sum = lhs + rhs;
    const auto diff = lhs - rhs;

    EXPECT_DOUBLE_EQ(sum[0][0], 1.5);
    EXPECT_DOUBLE_EQ(sum[1][1], 4.5);
    EXPECT_DOUBLE_EQ(diff[0][0], 0.5);
    EXPECT_DOUBLE_EQ(diff[1][1], 3.5);
}

TEST(CoreMatrix, SupportsScalarArithmeticAndConversion)
{
    uv::core::Matrix<double> m{1, 3};
    m[0][0] = 2.0;
    m[0][1] = 4.0;
    m[0][2] = 8.0;

    const auto scaled = 0.5 * (m * 4.0) / 2.0;
    const auto negated = -m;
    const auto asFloat = scaled.as<float>();

    EXPECT_DOUBLE_EQ(scaled[0][0], 2.0);
    EXPECT_DOUBLE_EQ(scaled[0][1], 4.0);
    EXPECT_DOUBLE_EQ(scaled[0][2], 8.0);
    EXPECT_DOUBLE_EQ(negated[0][2], -8.0);
    EXPECT_FLOAT_EQ(asFloat[0][1], 4.0F);
}
