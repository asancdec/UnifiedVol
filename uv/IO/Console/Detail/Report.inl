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

#include "Base/Macros/Inform.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"

#include <cstddef>
#include <format>
#include <iomanip>
#include <sstream>

namespace uv::io::report
{
namespace detail
{
template <typename HeaderVec, typename RowLabels, typename Matrix> void printMatrix(
    std::string_view title,
    const HeaderVec& header,
    const RowLabels& rowLabels,
    const Matrix& M,
    unsigned int headerPrec,
    unsigned int rowLabelPrec,
    unsigned int valuePrec
)
{
    std::ostringstream oss;
    oss << '\n' << title << '\t';

    oss << std::fixed << std::setprecision(precision(headerPrec));
    for (const auto& h : header)
        oss << h << '\t';
    oss << '\n';

    for (std::size_t i = 0; i < M.rows(); ++i)
    {
        oss << std::fixed << std::setprecision(precision(rowLabelPrec)) << rowLabels[i]
            << '\t';

        oss << std::fixed << std::setprecision(precision(valuePrec));
        for (const auto& v : M[i])
            oss << v << '\t';

        oss << '\n';
    }

    INFO(oss.str());
}

template <typename Vector> void printVector(const Vector& v, unsigned int valuePrec)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision(valuePrec));

    for (const auto& x : v)
        oss << x << '\t';

    oss << '\n';

    INFO(oss.str());
}
} // namespace detail

template <std::floating_point T>
void volatility(const core::VolSurface<T>& volSurface, unsigned int valuePrec)
{

    detail::printMatrix(
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
void volatility(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    volatility(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void totalVariance(const core::VolSurface<T>& volSurface, unsigned int valuePrec)
{

    detail::printMatrix(
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
void totalVariance(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    totalVariance(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void variance(const core::VolSurface<T>& volSurface, unsigned int valuePrec)
{

    detail::printMatrix(
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
void variance(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    variance(marketState.volSurface, valuePrec);
}

template <std::floating_point T>
void logKF(const core::VolSurface<T>& volSurface, unsigned int valuePrec)
{

    detail::printMatrix(
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
void logKF(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    logKF(marketState.volSurface, valuePrec);
}

template <std::floating_point T> void callPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec
)
{
    detail::printMatrix(
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
void callPrices(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    const core::VolSurface<T>& volSurface{marketState.volSurface};

    detail::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, marketState.interestCurve, true),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T> void putPrices(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    unsigned int valuePrec
)
{
    detail::printMatrix(
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
void putPrices(const core::MarketState<T>& marketState, unsigned int valuePrec)
{
    const core::VolSurface<T>& volSurface{marketState.volSurface};

    detail::printMatrix(
        "T\\K/S",
        volSurface.moneyness(),
        volSurface.maturities(),
        math::black::priceB76(volSurface, marketState.interestCurve, false),
        2,
        2,
        valuePrec
    );
}

template <std::floating_point T>
void sviParams(const models::svi::Params<T>& params, unsigned int valuePrec)
{
    INFO(std::format(
        "T={:.4f}, a={:.{}f}, b={:.{}f}, rho={:.{}f}, m={:.{}f}, sigma={:.{}f}",
        params.t,
        params.a,
        valuePrec,
        params.b,
        valuePrec,
        params.rho,
        valuePrec,
        params.m,
        valuePrec,
        params.sigma,
        valuePrec
    ));
}
} // namespace uv::io::report
