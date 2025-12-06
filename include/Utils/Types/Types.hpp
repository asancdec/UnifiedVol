/*
* Types.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <complex>
#include <vector>

namespace uv
{
    
    template <typename T>
    using Vector = std::vector<T>;

    template <typename T>
    using Matrix = std::vector<std::vector<T>>;

    using cplx = std::complex<long double>;

    template <typename From, typename To>
    Vector<To> convertVector(const Vector<From>& input);
}