// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::pde
{
template <std::floating_point T, std::size_t N> struct Grid
{
    std::array<T, N> x_;
    std::array<T, N - 1> dx_;
    std::array<T, N - 1> dx2_;

    Grid() = delete;

    explicit Grid(std::span<const T, N>);
};

} // namespace uv::math::pde

#include "Math/PDE/Detail/Grid.inl"
