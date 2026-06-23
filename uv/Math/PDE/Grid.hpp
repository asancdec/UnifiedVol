// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::pde
{

template <std::floating_point T, std::size_t N> class Grid
{
  private:
    std::array<T, N> x_;
    std::array<T, N - 1> dx_;

    void validate() const;
    void setGridSteps() noexcept;

  public:
    Grid() = delete;
    explicit Grid(std::span<const T, N>);

    std::span<const T> x() const noexcept;
    std::span<const T> dx() const noexcept;
};

template <std::floating_point T, std::size_t N>
Grid<T, N> generateCenteredSinHGrid(T xMin, T xMax, T beta = 0);

template <std::floating_point T, std::size_t N>
Grid<T, N> generateUniformGrid(T xMin, T xMax);

} // namespace uv::math::pde

#include "Math/PDE/Detail/Grid.inl"
