// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <span>

namespace uv::math::interp::bspline
{

template <std::floating_point T, int K = 3, bool doValidate = true> class BSpline
{
  public:
    BSpline() = delete;

    explicit BSpline(std::span<const T> cPoints, std::span<const T> knots);

    void evalInplace(std::span<T> out, std::span<const T> x) const;

    Vector<T> eval(std::span<const T> x) const noexcept;

  private:
    Vector<T> cPoints_;
    Vector<T> knots_;

    void validate() const;
};
} // namespace uv::math::interp::bspline

#include "Math/Interpolation/BSpline/Detail/Interpolator.inl"
