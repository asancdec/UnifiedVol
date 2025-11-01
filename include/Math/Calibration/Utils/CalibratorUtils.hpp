/**
* CalibratorUtils.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include <array>
#include <cstddef>
#include <string_view>

namespace uv
{   
    // Clamp initial guess within upper and lower bounds
    template <::std::size_t N>
    void clamp(::std::array<double, N>& initGuess,
        const ::std::array<double, N>& lowerBounds,
        const ::std::array<double, N>& upperBounds,
        const ::std::array<::std::string_view, N>& paramNames) noexcept;
}

#include "CalibratorUtils.inl"