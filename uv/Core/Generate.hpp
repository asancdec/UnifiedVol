// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <span>

#include "Base/Types.hpp"
#include "Core/Curve.hpp"
#include "Core/MarketData.hpp"
#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"

namespace uv::core
{

template <std::floating_point T> MarketState<T> generateMarketState(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Matrix<T>& vol
);

template <std::floating_point T> VolSurface<T>
generateVolSurface(const VolSurface<T>& volSurface, const core::Matrix<T>& vol);
} // namespace uv::core

#include "Core/Detail/Generate.inl"
