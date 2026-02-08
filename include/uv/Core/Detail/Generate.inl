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

#include <cstddef>

namespace uv::core
{

template <std::floating_point T>
MarketState<T> generateMarketState(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Matrix<T>& vol
)
{
    Curve<T> interestCurve{detail::generateInterestCurve<T>(marketData, maturities)};

    Curve<T> dividendCurve{detail::generateDividendCurve<T>(marketData, maturities)};

    VolSurface<T> volSurface{detail::generateVolSurface<
        T>(marketData, maturities, moneyness, interestCurve, dividendCurve, vol)};

    return MarketState<T>{
        .interestCurve = std::move(interestCurve),
        .dividendCurve = std::move(dividendCurve),
        .volSurface = std::move(volSurface)
    };
}

template <std::floating_point T>
VolSurface<T>
generateVolSurface(const VolSurface<T>& volSurface, const core::Matrix<T>& vol)
{
    return VolSurface<T>{
        volSurface.maturities(),
        volSurface.forwards(),
        volSurface.strikes(),
        volSurface.moneyness(),
        vol
    };
}

} // namespace uv::core

namespace uv::core::detail
{
template <std::floating_point T>
Curve<T>
generateInterestCurve(const MarketData<T>& marketData, std::span<const T> maturities)
{
    return Curve{marketData.interestRate, maturities};
}

template <std::floating_point T>
Curve<T>
generateDividendCurve(const MarketData<T>& marketData, std::span<const T> maturities)
{
    return Curve{marketData.dividendYield, maturities};
}

template <std::floating_point T>
VolSurface<T> generateVolSurface(
    const MarketData<T>& marketData,
    std::span<const T> maturities,
    std::span<const T> moneyness,
    const Curve<T>& interestCurve,
    const Curve<T>& dividendCurve,
    const Matrix<T>& vol
)
{
    return VolSurface<T>{
        maturities,
        generateForwards(marketData.spot, maturities, interestCurve, dividendCurve),
        generateStrikes(marketData.spot, moneyness),
        moneyness,
        vol

    };
}

template <std::floating_point T>
Vector<T> generateForwards(
    const T spot,
    std::span<const T> maturities,
    const Curve<T>& interestCurve,
    const Curve<T>& dividendCurve
)
{
    const std::size_t n{maturities.size()};

    Vector<T> forwards;
    forwards.resize(n);

    for (std::size_t i{0}; i < n; ++i)
    {
        const T maturity{maturities[i]};
        forwards[i] = spot * dividendCurve.interpolateDF(maturity) /
                      interestCurve.interpolateDF(maturity);
    }

    return forwards;
}

template <std::floating_point T>
Vector<T> generateStrikes(T spot, const std::span<const T> moneyness) noexcept
{
    const std::size_t n{moneyness.size()};

    Vector<T> strikes;
    strikes.resize(n);

    for (std::size_t i{0}; i < n; ++i)
    {
        strikes[i] = moneyness[i] * spot;
    }

    return strikes;
}
} // namespace uv::core::detail
