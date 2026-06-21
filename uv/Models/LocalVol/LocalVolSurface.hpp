// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Math/Interpolation/BSpline/Interpolator.hpp"

#include <concepts>
#include <span>

namespace uv::models::local_vol
{

template <std::floating_point T> class LocalVolSlice
{
  public:
    LocalVolSlice() = delete;

    LocalVolSlice(std::span<const T> cPoints, std::span<const T> knots);

    void localVarInplace(std::span<T>, std::span<const T>) const;
    Vector<T> localVar(std::span<const T>) const;

    void localVolInplace(std::span<T>, std::span<const T>) const;
    Vector<T> localVol(std::span<const T>) const;

  private:
    uv::math::interp::bspline::BSpline<T> interpolator_;

    Vector<T> exponentiateControlPoints(std::span<const T>) const;
};

} // namespace uv::models::local_vol

#include "Models/LocalVol/Detail/LocalVolSurface.inl"
