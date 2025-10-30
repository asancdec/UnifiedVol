/**
* MathFunctions.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include <complex>

namespace uv  
{
    /**
     * @brief Numerically stable evaluation of log(1 + z) for complex arguments.
     *
     * This function avoids cancellation when |z| is small by switching to a
     * truncated Taylor expansion. The full definition is provided in
     * MathFunctions.inl.
     *
     * @tparam T Floating-point type (float, double, long double)
     * @param z Complex argument
     * @return std::complex<T> Accurate value of log(1 + z)
     */
    template <typename T>
    std::complex<T> log1pComplex(const std::complex<T>& z) noexcept;
}

#include "MathFunctions.inl"