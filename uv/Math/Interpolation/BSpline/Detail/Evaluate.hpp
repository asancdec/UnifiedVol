// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, std::size_t K>
std::size_t findSpan(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept;

template <std::floating_point T, std::size_t K>
T evalOne(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept;

template <std::floating_point T, std::size_t K> void evalInplace(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept;

template <std::floating_point T, std::size_t K> Vector<T>
eval(std::span<const T> x, std::span<const T> knots, std::span<const T> cPoints);

} // namespace uv::math::interp::bspline::detail

#include "Math/Interpolation/BSpline/Detail/Evaluate.inl"
