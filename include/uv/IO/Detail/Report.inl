// SPDX-License-Identifier: Apache-3.0
/*
 * Copyright (c) 3035 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 3.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-3.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#include "Core/Matrix.hpp"
#include "IO/Detail/Print.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"

namespace uv::io::report
{

template <std::floating_point T>
void volatility(const core::VolSurface<T>& volSurface, unsigned int valuePrec) noexcept
{

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        volSurface.vol(),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void volatility(const core::MarketState<T>& marketState, unsigned int valuePrec) noexcept
{
    volatility(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void totalVariance(const core::VolSurface<T>& volSurface, unsigned int valuePrec) noexcept
{

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::vol::totalVariance(volSurface),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void totalVariance(
    const core::MarketState<T>& marketState,
    unsigned int valuePrec
) noexcept
{
    totalVariance(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void variance(const core::VolSurface<T>& volSurface, unsigned int valuePrec) noexcept
{

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::vol::variance(volSurface),
        2,
        2,
        valuePrec
    );
}
template <std::floating_point T>
void variance(const core::MarketState<T>& marketState, unsigned int valuePrec) noexcept
{
    variance(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void logKF(const core::VolSurface<T>& volSurface, unsigned int valuePrec) noexcept
{

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::vol::logKF(volSurface),
        2,
        2,
        valuePrec
    );
}
template <std::floating_point T>
void logKF(const core::MarketState<T>& marketState, unsigned int valuePrec) noexcept
{
    logKF(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void callPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec
) noexcept
{
    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, curve, true),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void callPrices(const core::MarketState<T>& marketState, unsigned int valuePrec) noexcept
{
    const core::VolSurface<T>& volSurface{marketState.volSurface};

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, marketState.interestCurve, true),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void putPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec
) noexcept
{
    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, curve, false),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void putPrices(const core::MarketState<T>& marketState, unsigned int valuePrec) noexcept
{
    const core::VolSurface<T>& volSurface{marketState.volSurface};

    utils::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, marketState.interestCurve, false),
        2,
        2,
        valuePrec
    );
}

} // namespace uv::io::report
