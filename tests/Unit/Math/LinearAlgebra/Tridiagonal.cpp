// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/Tridiagonal.hpp"
#include "Support/Math/LinearAlgebra/Tridiagonal.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <random>

namespace uv::tests::unit::linear_algebra::detail
{
template <typename T, std::size_t N> void solveAndExpectNear(
    std::array<T, N> rhs,
    const std::array<T, N>& upper,
    const std::array<T, N>& middle,
    const std::array<T, N>& lower,
    const std::array<T, N>& expected,
    const double tolerance
)
{
    std::array<T, N> scratch{};

    uv::math::linear_algebra::thomasSolve<T, N>(rhs, upper, middle, lower, scratch);

    for (std::size_t i{0}; i < N; ++i)
    {
        EXPECT_NEAR(rhs[i], expected[i], tolerance) << "i=" << i;
    }
}

template <std::size_t N> void runRandomDiagonallyDominantCases(const std::size_t cases)
{
    std::mt19937_64 rng{0x5A17D1A6D0ULL + N};
    std::uniform_real_distribution<double> offDiagonal{-0.75, 0.75};
    std::uniform_real_distribution<double> solution{-5.0, 5.0};
    std::uniform_real_distribution<double> margin{0.25, 2.0};

    for (std::size_t testCase{0}; testCase < cases; ++testCase)
    {
        std::array<double, N> upper{};
        std::array<double, N> middle{};
        std::array<double, N> lower{};
        std::array<double, N> expected{};

        for (std::size_t i{0}; i < N; ++i)
        {
            upper[i] = (i + 1 < N) ? offDiagonal(rng) : 0.0;
            lower[i] = (i > 0) ? offDiagonal(rng) : 0.0;
            middle[i] = std::abs(upper[i]) + std::abs(lower[i]) + margin(rng);
            expected[i] = solution(rng);
        }

        const auto rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
            upper,
            middle,
            lower,
            expected
        )};

        solveAndExpectNear(rhs, upper, middle, lower, expected, 1e-11);
    }
}
} // namespace uv::tests::unit::linear_algebra::detail

using namespace uv::tests::unit::linear_algebra::detail;

TEST(MathTridiagonal, ThomasSolveRecoversKnownSolution)
{
    constexpr std::size_t n{5};
    const std::array<double, n> upper{-1.0, -1.0, -1.0, -1.0, 0.0};
    const std::array<double, n> middle{4.0, 4.0, 4.0, 4.0, 4.0};
    const std::array<double, n> lower{0.0, -1.0, -1.0, -1.0, -1.0};
    const std::array<double, n> expected{1.0, 2.0, 3.0, 4.0, 5.0};
    const std::array<double, n> rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
        upper,
        middle,
        lower,
        expected
    )};

    solveAndExpectNear(rhs, upper, middle, lower, expected, 1e-14);
}

TEST(MathTridiagonal, ThomasSolveHandlesTwoByTwoSystem)
{
    constexpr std::size_t n{2};
    const std::array<double, n> upper{2.0, 0.0};
    const std::array<double, n> middle{5.0, 7.0};
    const std::array<double, n> lower{0.0, 3.0};
    const std::array<double, n> expected{2.0, -1.0};
    const std::array<double, n> rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
        upper,
        middle,
        lower,
        expected
    )};

    solveAndExpectNear(rhs, upper, middle, lower, expected, 1e-14);
}

TEST(MathTridiagonal, ThomasSolveHandlesRandomDiagonallyDominantSystems)
{
    runRandomDiagonallyDominantCases<3>(100);
    runRandomDiagonallyDominantCases<16>(100);
    runRandomDiagonallyDominantCases<64>(50);
}

TEST(MathTridiagonal, ThomasSolveIgnoresUnusedEndpointCoefficients)
{
    constexpr std::size_t n{4};
    const std::array<double, n> upper{0.4, -0.3, 0.2, 1.0e12};
    const std::array<double, n> middle{2.5, 3.0, 3.5, 4.0};
    const std::array<double, n> lower{-1.0e12, -0.2, 0.5, -0.4};
    const std::array<double, n> expected{1.0, -2.0, 0.5, 3.0};
    const std::array<double, n> rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
        upper,
        middle,
        lower,
        expected
    )};

    solveAndExpectNear(rhs, upper, middle, lower, expected, 1e-13);
}

TEST(MathTridiagonal, ThomasSolveSupportsFloat)
{
    constexpr std::size_t n{4};
    const std::array<float, n> upper{-0.5F, 0.25F, -0.1F, 0.0F};
    const std::array<float, n> middle{2.0F, 2.5F, 3.0F, 3.5F};
    const std::array<float, n> lower{0.0F, 0.4F, -0.2F, 0.3F};
    const std::array<float, n> expected{1.25F, -0.5F, 2.0F, 0.75F};
    const std::array<float, n> rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
        upper,
        middle,
        lower,
        expected
    )};

    solveAndExpectNear(rhs, upper, middle, lower, expected, 1e-5);
}
