// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::linear_algebra
{
template <std::floating_point T, std::size_t N> void thomasSolve(
    std::span<T, N> x,
    std::span<const T, N> upper,
    std::span<const T, N> middle,
    std::span<const T, N> lower,
    std::span<T, N> scratch
) noexcept;
} // namespace uv::math::linear_algebra

#include "Math/LinearAlgebra/Detail/Tridiagonal.inl"
