// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/VolSurface.hpp"

#include <concepts>

namespace uv::core
{
template <std::floating_point T> struct MarketState
{
    Curve<T> interestCurve;
    Curve<T> dividendCurve;
    VolSurface<T> volSurface;
};

} // namespace uv::core