// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/MatrixOps.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixOps, GenerateAndTransformIndexed)
{
    const auto m = uv::math::linear_algebra::generateIndexed<double>(
        2,
        3,
        [](std::size_t i, std::size_t j)
        {
            return static_cast<double>(10 * i + j);
        }
    );

    const auto transformed = uv::math::linear_algebra::transformIndexed<double>(
        m,
        [](std::size_t i, std::size_t j, double x)
        {
            return x + static_cast<double>(i + j);
        }
    );

    EXPECT_DOUBLE_EQ(m[1][2], 12.0);
    EXPECT_DOUBLE_EQ(transformed[1][2], 15.0);
}

TEST(MathMatrixOps, HadamardDivideSquareSqrtAndReciprocal)
{
    uv::core::Matrix<double> lhs{2, 2};
    lhs[0][0] = 2.0;
    lhs[0][1] = 4.0;
    lhs[1][0] = 8.0;
    lhs[1][1] = 16.0;

    uv::core::Matrix<double> rhs{2, 2, 2.0};

    const auto product = uv::math::linear_algebra::hadamard(lhs, rhs);
    const auto quotient = uv::math::linear_algebra::divide(product, rhs);
    const auto squared = uv::math::linear_algebra::square(rhs);
    const auto roots = uv::math::linear_algebra::sqrt(squared);
    const auto reciprocal = uv::math::linear_algebra::reciprocal(rhs);

    EXPECT_DOUBLE_EQ(product[1][1], 32.0);
    EXPECT_DOUBLE_EQ(quotient[1][1], lhs[1][1]);
    EXPECT_DOUBLE_EQ(squared[0][0], 4.0);
    EXPECT_DOUBLE_EQ(roots[0][0], 2.0);
    EXPECT_DOUBLE_EQ(reciprocal[0][0], 0.5);
}

TEST(MathMatrixOps, HadamardWithVectorScalesRows)
{
    uv::core::Matrix<double> m{2, 2, 3.0};
    const std::vector<double> rowScale{2.0, 4.0};

    const auto scaled = uv::math::linear_algebra::hadamard(m, rowScale);

    EXPECT_DOUBLE_EQ(scaled[0][0], 6.0);
    EXPECT_DOUBLE_EQ(scaled[0][1], 6.0);
    EXPECT_DOUBLE_EQ(scaled[1][0], 12.0);
    EXPECT_DOUBLE_EQ(scaled[1][1], 12.0);
}
