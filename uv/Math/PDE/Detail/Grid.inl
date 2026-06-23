// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"

#include <algorithm>
#include <cmath>

namespace uv::math::pde::detail
{
template <std::floating_point T, std::size_t N>
std::array<T, N> generateUniformGrid(T xMin, T xMax);
} // namespace uv::math::pde::detail

namespace uv::math::pde
{
template <std::floating_point T, std::size_t N> Grid<T, N>::Grid(std::span<const T, N> x)
{
    std::copy(x.begin(), x.end(), x_.begin());
    validate();
    setGridSteps();
}

template <std::floating_point T, std::size_t N> void Grid<T, N>::validate() const
{
    static_assert(N > 2, "Size of grid must be larger than 2");
    REQUIRE_STRICTLY_MONOTONIC(x_);
}

template <std::floating_point T, std::size_t N> void Grid<T, N>::setGridSteps() noexcept
{
    for (std::size_t i{0}; i < N - 1; ++i)
    {
        dx_[i] = x_[i + 1] - x_[i];
    }
}

template <std::floating_point T, std::size_t N>
std::span<const T> Grid<T, N>::x() const noexcept
{
    return x_;
}

template <std::floating_point T, std::size_t N>
std::span<const T> Grid<T, N>::dx() const noexcept
{
    return dx_;
}

template <std::floating_point T, std::size_t N>
Grid<T, N> generateCenteredSinHGrid(T xMin, T xMax, T beta)
{
    constexpr T uniformThreshold{1e-10};

    REQUIRE_NON_NEGATIVE(beta);
    REQUIRE_CLOSE(xMin, -xMax, T{1e-12});

    if (beta < std::abs(uniformThreshold))
        return generateUniformGrid<T, N>(xMin, xMax);

    std::array<T, N> x{detail::generateUniformGrid<T, N>(T{-1.0}, T{1.0})};

    const T scale{xMax / std::sinh(beta)};

    for (std::size_t i{0}; i < N; ++i)
    {
        x[i] = scale * std::sinh(x[i] * beta);
    }

    return Grid<T, N>{x};
}

template <std::floating_point T, std::size_t N>
Grid<T, N> generateUniformGrid(T xMin, T xMax)
{
    return Grid<T, N>{detail::generateUniformGrid<T, N>(xMin, xMax)};
} // namespace uv::math::pde::detail

} // namespace uv::math::pde

namespace uv::math::pde::detail
{
template <std::floating_point T, std::size_t N>
std::array<T, N> generateUniformGrid(T xMin, T xMax)
{
    REQUIRE_GREATER(xMax, xMin);
    static_assert(N > 2, "Size of grid must be larger than 2");

    const T range{xMax - xMin};
    const T stepSize{range / static_cast<T>(N - 1)};

    std::array<T, N> x{};

    for (std::size_t i{0}; i < N; ++i)
    {
        x[i] = xMin + stepSize * static_cast<T>(i);
    }

    return x;
}
} // namespace uv::math::pde::detail
