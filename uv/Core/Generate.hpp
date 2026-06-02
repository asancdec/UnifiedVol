// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <span>
#include <utility>

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

namespace detail
{

template <std::floating_point T> Curve<T>
generateInterestCurve(const MarketData<T>& marketData, std::span<const T> maturities);

template <std::floating_point T> Curve<T>
generateDividendCurve(const MarketData<T>& marketData, std::span<const T> maturities);

template <std::floating_point T> VolSurface<T> generateVolSurface(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Curve<T>& interestCurve,
    const Curve<T>& dividendCurve,
    const Matrix<T>& vol
);

template <std::floating_point T> Vector<T> generateForwards(
    const T spot,
    std::span<const T> maturities,
    const Curve<T>& interestCurve,
    const Curve<T>& dividendCurve
);

template <std::floating_point T>
Vector<T> generateStrikes(T spot, const std::span<const T> moneyness) noexcept;

} // namespace detail
} // namespace uv::core

#include "Core/Detail/Generate.inl"
