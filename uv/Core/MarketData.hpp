// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>

namespace uv::core
{

template <std::floating_point T> struct MarketData
{
    T interestRate;
    T dividendYield;
    T spot;
};
} // namespace uv::core