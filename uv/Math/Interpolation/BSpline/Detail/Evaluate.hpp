// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, int K> Vector<T>
eval(std::span<const T> x, std::span<const T> knots, std::span<const T> cPoints) noexcept;

template <std::floating_point T, int K> void evalInplace(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept;

template <std::floating_point T, int K> void coxDeBoor(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::size_t idx
) noexcept;

template <std::floating_point T> void basisZero(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::size_t idx
) noexcept;

} // namespace uv::math::interp::bspline::detail

#include "Math/Interpolation/BSpline/Detail/Evaluate.inl"
