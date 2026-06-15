// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/Tridiagonal.hpp"
#include "Support/Math/LinearAlgebra/Tridiagonal.hpp"
#include "Support/Performance/Budgets.hpp"
#include "Support/Performance/Timing.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>

namespace
{
template <std::size_t N> std::array<double, N> makeExpectedSolution()
{
    std::array<double, N> x{};
    for (std::size_t i{0}; i < N; ++i)
    {
        const double t{static_cast<double>(i)};
        x[i] = 1.0 + 0.01 * std::sin(0.07 * t) + 0.02 * std::cos(0.03 * t);
    }
    return x;
}

} // namespace

TEST(PerformanceTridiagonal, SolvesRepeatedSystemsWithinLatencyBudget)
{
    const auto budget = uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        uv::tests::performance::TridiagonalThomasSolveBudgetKey
    );

    constexpr std::size_t n{512};
    constexpr std::size_t iterations{2'000};

    std::array<double, n> upper{};
    std::array<double, n> middle{};
    std::array<double, n> lower{};

    for (std::size_t i{0}; i < n; ++i)
    {
        upper[i] = (i + 1 < n) ? -0.5 : 0.0;
        middle[i] = 2.0;
        lower[i] = (i > 0) ? -0.5 : 0.0;
    }

    const auto expected{makeExpectedSolution<n>()};
    const auto rhs{uv::tests::math::linear_algebra::multiplyTridiagonal(
        upper,
        middle,
        lower,
        expected
    )};

    std::array<double, n> x{};
    std::array<double, n> scratch{};
    double checksum{};

    const double ms = uv::tests::performance::bestElapsedMs(
        [&]
        {
            checksum = 0.0;
            for (std::size_t i{0}; i < iterations; ++i)
            {
                x = rhs;
                uv::math::linear_algebra::thomasSolve<double, n>(
                    x,
                    upper,
                    middle,
                    lower,
                    scratch
                );
                checksum += x[i % n];
            }
        }
    );

    EXPECT_NEAR(x[n / 2], expected[n / 2], 1e-12);
    EXPECT_TRUE(std::isfinite(checksum));
    EXPECT_LT(ms, budget.maxMs);
}
