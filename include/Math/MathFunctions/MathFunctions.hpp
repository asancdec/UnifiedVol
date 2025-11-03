/**
* MathFunctions.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <complex>
#include <limits>

namespace uv  
{
    // Numerically stable evaluation of log(1 + z) for complex arguments.
    template <typename T>
    ::std::complex<T> log1pComplex(const ::std::complex<T>& z) noexcept;

    // cosm1(b) := cos(b) - 1 = -2 sin²(b / 2)
    template <typename T>
    T cosm1(T b) noexcept;

    // expm1Complex(z) := e^z - 1 (stable for small |z|)
    template <typename T>
    ::std::complex<T> expm1Complex(const ::std::complex<T>& z) noexcept;

    // Normal cumulative density function
    template <typename T>
    T normalCDF(T x) noexcept;

    // Normal probability density function
    template <typename T>
    T normalPDF(T x) noexcept;

    // Black-Scholes pricing
    template <typename T>
    T blackScholes(T t,
        T r,
        T q,
        T vol,
        T S,
        T K,
        bool isCall = true) noexcept;

    // Black-Scholes Vega
    template <typename T>
    T vegaBS(T d1,
        T t,
        T q,
        T S) noexcept;

    // Black-Scholes Volga
    template <typename T>
    T volgaBS(T vega,
        T d1,
        T t,
        T vol) noexcept;

    // Calculate implied volatility using Halley's method
    double impliedVolBS(double mktPriceBS,
        double t,
        double r,
        double q,
        double S,
        double K,
        bool isCall = true,
        double ftolAbs = ::std::numeric_limits<double>::epsilon(),
        unsigned int maxEval = 100) noexcept;
}

#include "MathFunctions.inl"