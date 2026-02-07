// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#pragma once

#include <concepts>
#include <span>
#include <string>
#include <utility>

#include "Base/Alias.hpp"
#include "Core/Curve.hpp"
#include "Core/MarketData.hpp"
#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"

namespace uv::core
{

template <std::floating_point T>
MarketState<T> generateMarketState(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Matrix<T>& vol
);

namespace detail
{

template <std::floating_point T>
Curve<T>
generateInterestCurve(const MarketData<T>& marketData, std::span<const T> maturities);

template <std::floating_point T>
Curve<T>
generateDividendCurve(const MarketData<T>& marketData, std::span<const T> maturities);

template <std::floating_point T>
VolSurface<T> generateVolSurface(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Curve<T>& interestCurve,
    const Curve<T>& dividendCurve,
    const Matrix<T>& vol
);

template <std::floating_point T>
Vector<T> generateForwards(
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