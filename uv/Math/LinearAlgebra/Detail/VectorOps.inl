// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Base/Types.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <numeric>
#include <span>
#include <utility>
#include <vector>

namespace uv::math::linear_algebra
{
template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept
{
    evalInplace<T, N, F>(grid, std::forward<F>(f));
    return grid;
}

template <std::floating_point T, std::size_t N, typename F>
void evalInplace(std::array<T, N>& grid, F&& f) noexcept
{
    auto&& func = std::forward<F>(f);

    for (std::size_t i{0}; i < N; ++i)
        grid[i] = std::invoke(func, grid[i]);
}

template <std::floating_point T> T sum(std::span<const T> x) noexcept
{
    return std::accumulate(x.begin(), x.end(), T{});
}

template <std::floating_point T>
Vector<T> multiply(std::span<const T> v, const T x) noexcept
{
    std::size_t vSize{v.size()};

    Vector<T> result(vSize);

    for (std::size_t i = 0; i < vSize; ++i)
    {
        result[i] = v[i] * x;
    }

    return result;
}

template <std::floating_point T> Vector<T> reciprocal(std::span<const T> v) noexcept
{
    std::size_t vSize{v.size()};

    Vector<T> result(vSize);

    for (std::size_t i = 0; i < vSize; ++i)
    {
        result[i] = T(1.0) / v[i];
    }

    return result;
}

template <std::floating_point T>
Vector<T> hadamard(std::span<const T> a, std::span<const T> b)
{

    REQUIRE_SAME_SIZE(a, b);

    std::size_t aSize{a.size()};

    Vector<T> c(aSize);

    for (std::size_t i = 0; i < aSize; ++i)
    {
        c[i] = a[i] * b[i];
    }

    return c;
}

template <typename T> Vector<T> makeSequence(std::size_t n, T start) noexcept
{
    Vector<T> v(n);
    std::iota(v.begin(), v.end(), start);
    return v;
}

template <typename T> T minValue(std::span<const T> x)
{
    REQUIRE_NON_EMPTY(x);

    return *std::min_element(x.begin(), x.end());
}

template <typename T> T maxValue(std::span<const T> x)
{
    REQUIRE_NON_EMPTY(x);

    return *std::max_element(x.begin(), x.end());
}

} // namespace uv::math::linear_algebra
