// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::interp::bspline
{

template <std::floating_point T, std::size_t K = 3, bool doValidate = true> class BSpline
{
  public:
    BSpline() = delete;

    explicit BSpline(std::span<const T> cPoints, std::span<const T> knots);

    void evalInplace(std::span<T>, std::span<const T>) const;

    Vector<T> eval(std::span<const T>) const;

    void setControlPoints(std::span<const T>);

    void setKnots(std::span<const T>);

  private:
    Vector<T> cPoints_;
    Vector<T> knots_;
};
} // namespace uv::math::interp::bspline

#include "Math/Interpolation/BSpline/Detail/Interpolator.inl"
