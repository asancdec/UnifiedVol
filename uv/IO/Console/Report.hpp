// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/Params.hpp"

#include <concepts>
#include <string_view>

namespace uv::io::report
{
namespace detail
{
int precision(unsigned int) noexcept;

template <typename HeaderVec, typename RowLabels, typename Matrix> void printMatrix(
    std::string_view title,
    const HeaderVec& header,
    const RowLabels& rowLabels,
    const Matrix& M,
    unsigned int headerPrec = 2,
    unsigned int rowLabelPrec = 2,
    unsigned int valuePrec = 5
);

template <typename Vector> void printVector(const Vector& v, unsigned int valuePrec = 5);
} // namespace detail

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

#include "IO/Console/Detail/Report.inl"
