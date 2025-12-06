/*
* Utils.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <algorithm>
#include <string> 
#include <array>
#include <cstddef>
#include <string_view>
#include <span>

namespace uv::cal
{
    // Clamp initial guess within upper and lower bounds
    template <std::size_t N>
    void clamp(std::array<Real, N>& initGuess,
        const std::array<Real, N>& lowerBounds,
        const std::array<Real, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Warn if upper or lower bounds are touched
    template <std::size_t N>
    void warnBoundsHit(std::span<Real> x,
        const std::array<Real, N>& lowerBounds,
        const std::array<Real, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Log calibration results 
    template <std::size_t N>
    void logResults(std::span<Real> x,
        const std::array<std::string_view, N>& paramNames,
        Real sse,
        unsigned iterCount,
        Real elapsedMs,
        bool isSuccess) noexcept;
}

#include "Utils.inl"