// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Inform.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"

#include <algorithm>
#include <cstddef>
#include <format>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>

namespace uv::io::report
{
namespace detail
{
[[nodiscard]] inline int precision(unsigned int value) noexcept
{
    const unsigned int capped{
        std::min(value, static_cast<unsigned int>(std::numeric_limits<int>::max()))
    };
    return static_cast<int>(capped);
}

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
    std::string output;
    std::format_to(std::back_inserter(output), "\n{}\t", title);

    for (const auto& h : header)
        std::format_to(std::back_inserter(output), "{:.{}f}\t", h, precision(headerPrec));
    output.push_back('\n');

    for (std::size_t i = 0; i < M.rows(); ++i)
    {
        std::format_to(
            std::back_inserter(output),
            "{:.{}f}\t",
            rowLabels[i],
            precision(rowLabelPrec)
        );

        for (const auto& v : M[i])
            std::format_to(
                std::back_inserter(output),
                "{:.{}f}\t",
                v,
                precision(valuePrec)
            );

        output.push_back('\n');
    }

    INFO(output);
}

template <typename Vector> void printVector(const Vector& v, unsigned int valuePrec)
{
    std::string output;

    for (const auto& x : v)
        std::format_to(std::back_inserter(output), "{:.{}f}\t", x, precision(valuePrec));

    output.push_back('\n');

    INFO(output);
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
