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

#include <Core/MarketState.hpp>
#include <Core/VolSurface.hpp>

#include <concepts>

namespace uv::io::report
{

template <std::floating_point T>
void volatility(
    const core::VolSurface<T>& volSurface,
    unsigned int valuePrec = 5
) noexcept;

template <std::floating_point T>
void volatility(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec = 5
) noexcept;

template <std::floating_point T>
void totalVariance(
    const core::VolSurface<T>& volSurface,
    unsigned int valuePrec = 5
) noexcept;

template <std::floating_point T>
void totalVariance(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec = 5
) noexcept;

template <std::floating_point T>
void variance(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 5) noexcept;

template <std::floating_point T>
void variance(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec = 5
) noexcept;

template <std::floating_point T>
void logKF(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 4) noexcept;

template <std::floating_point T>
void logKF(const core::MarketState<T>& marketState, unsigned int valuePrec = 4) noexcept;

template <std::floating_point T>
void callPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec = 3
) noexcept;

template <std::floating_point T>
void callPrices(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec = 3
) noexcept;

template <std::floating_point T>
void putPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec = 3
) noexcept;

template <std::floating_point T>
void putPrices(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec = 3
) noexcept;

} // namespace uv::io::report

#include <IO/Detail/Report.inl>