// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/VectorOps.hpp"

#include <array>
#include <gtest/gtest.h>
#include <vector>

TEST(MathVectorOps, EvaluatesAndReducesVectors)
{
    std::array<double, 3> grid{1.0, 2.0, 3.0};

    const auto squared = uv::math::linear_algebra::eval<double, 3>(
        grid,
        [](double x)
        {
            return x * x;
        }
    );

    EXPECT_DOUBLE_EQ(squared[0], 1.0);
    EXPECT_DOUBLE_EQ(squared[1], 4.0);
    EXPECT_DOUBLE_EQ(squared[2], 9.0);
    EXPECT_DOUBLE_EQ(uv::math::linear_algebra::sum<double>(squared), 14.0);
}

TEST(MathVectorOps, ElementwiseHelpersReturnExpectedValues)
{
    const std::vector<double> a{2.0, 4.0, 8.0};
    const std::vector<double> b{0.5, 0.25, 0.125};

    const auto scaled = uv::math::linear_algebra::multiply<double>(a, 2.0);
    const auto inv = uv::math::linear_algebra::reciprocal<double>(a);
    const auto product = uv::math::linear_algebra::hadamard<double>(a, b);

    EXPECT_DOUBLE_EQ(scaled[2], 16.0);
    EXPECT_DOUBLE_EQ(inv[0], 0.5);
    EXPECT_DOUBLE_EQ(product[0], 1.0);
    EXPECT_DOUBLE_EQ(product[2], 1.0);
}

TEST(MathVectorOps, SequenceAndExtremaHelpersWork)
{
    const auto sequence = uv::math::linear_algebra::makeSequence<int>(4, 3);

    EXPECT_EQ(sequence, (std::vector<int>{3, 4, 5, 6}));
    EXPECT_EQ(uv::math::linear_algebra::minValue<int>(sequence), 3);
    EXPECT_EQ(uv::math::linear_algebra::maxValue<int>(sequence), 6);
}
