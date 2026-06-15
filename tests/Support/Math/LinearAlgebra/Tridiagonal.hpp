// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstddef>

namespace uv::tests::math::linear_algebra
{
template <typename T, std::size_t N> std::array<T, N> multiplyTridiagonal(
    const std::array<T, N>& upper,
    const std::array<T, N>& middle,
    const std::array<T, N>& lower,
    const std::array<T, N>& x
)
{
    std::array<T, N> rhs{};

    for (std::size_t i{0}; i < N; ++i)
    {
        rhs[i] = middle[i] * x[i];
        if (i > 0)
        {
            rhs[i] += lower[i] * x[i - 1];
        }
        if (i + 1 < N)
        {
            rhs[i] += upper[i] * x[i + 1];
        }
    }

    return rhs;
}
} // namespace uv::tests::math::linear_algebra
