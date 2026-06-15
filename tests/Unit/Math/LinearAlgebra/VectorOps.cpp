// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/VectorOps.hpp"
#include "Base/Errors/Errors.hpp"

#include <array>
#include <cmath>
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
    const auto exp =
        uv::math::linear_algebra::exponential<double>(std::vector<double>{0.0, 1.0});
    const auto product = uv::math::linear_algebra::hadamard<double>(a, b);

    EXPECT_DOUBLE_EQ(scaled[2], 16.0);
    EXPECT_DOUBLE_EQ(inv[0], 0.5);
    EXPECT_DOUBLE_EQ(exp[0], 1.0);
    EXPECT_DOUBLE_EQ(exp[1], std::exp(1.0));
    EXPECT_DOUBLE_EQ(product[0], 1.0);
    EXPECT_DOUBLE_EQ(product[2], 1.0);
}

TEST(MathVectorOps, HadamardRequiresSameSize)
{
    const std::vector<double> a{1.0, 2.0, 3.0};
    const std::vector<double> b{1.0, 2.0};

    EXPECT_THROW(
        uv::math::linear_algebra::hadamard<double>(a, b),
        uv::errors::UnifiedVolError
    );
}

TEST(MathVectorOps, SquareRootReturnsElementwiseRoots)
{
    const std::vector<double> x{0.0, 1.0, 4.0, 9.0, 16.0};

    const auto roots = uv::math::linear_algebra::squareRoot<double>(x);

    ASSERT_EQ(roots.size(), x.size());
    EXPECT_DOUBLE_EQ(roots[0], 0.0);
    EXPECT_DOUBLE_EQ(roots[1], 1.0);
    EXPECT_DOUBLE_EQ(roots[2], 2.0);
    EXPECT_DOUBLE_EQ(roots[3], 3.0);
    EXPECT_DOUBLE_EQ(roots[4], 4.0);
}

TEST(MathVectorOps, SquareRootInplaceWritesOutput)
{
    const std::vector<double> x{2.25, 6.25, 12.25};
    std::vector<double> out(x.size());

    uv::math::linear_algebra::squareRootInplace<double>(out, x);

    EXPECT_DOUBLE_EQ(out[0], 1.5);
    EXPECT_DOUBLE_EQ(out[1], 2.5);
    EXPECT_DOUBLE_EQ(out[2], 3.5);
}

TEST(MathVectorOps, SquareRootInplaceRequiresSameSize)
{
    const std::vector<double> x{1.0, 4.0, 9.0};
    std::vector<double> out(2);

    EXPECT_THROW(
        uv::math::linear_algebra::squareRootInplace<double>(out, x),
        uv::errors::UnifiedVolError
    );
}

TEST(MathVectorOps, SequenceAndExtremaHelpersWork)
{
    const auto sequence = uv::math::linear_algebra::makeSequence<int>(4, 3);

    EXPECT_EQ(sequence, (std::vector<int>{3, 4, 5, 6}));
    EXPECT_EQ(uv::math::linear_algebra::minValue<int>(sequence), 3);
    EXPECT_EQ(uv::math::linear_algebra::maxValue<int>(sequence), 6);
}
