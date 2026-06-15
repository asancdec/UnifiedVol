// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/Params.hpp"

#include <concepts>
#include <ranges>

namespace uv::io::report
{

template <std::floating_point T>
void volatility(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 5);

template <std::floating_point T>
void volatility(const core::MarketState<T>& marketState, unsigned int valuePrec = 5);

template <std::floating_point T>
void totalVariance(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 5);

template <std::floating_point T>
void totalVariance(const core::MarketState<T>& marketState, unsigned int valuePrec = 5);

template <std::floating_point T>
void variance(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 5);

template <std::floating_point T>
void variance(const core::MarketState<T>& marketState, unsigned int valuePrec = 5);

template <std::floating_point T>
void logKF(const core::VolSurface<T>& volSurface, unsigned int valuePrec = 4);

template <std::floating_point T>
void logKF(const core::MarketState<T>& marketState, unsigned int valuePrec = 4);

template <std::floating_point T> void callPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec = 3
);

template <std::floating_point T>
void callPrices(const core::MarketState<T>& marketState, unsigned int valuePrec = 3);

template <std::floating_point T> void putPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec = 3
);

template <std::floating_point T>
void putPrices(const core::MarketState<T>& marketState, unsigned int valuePrec = 3);

template <std::floating_point T>
void sviParams(const models::svi::Params<T>& params, unsigned int valuePrec = 6);

} // namespace uv::io::report

#include "IO/Detail/Report.inl"
