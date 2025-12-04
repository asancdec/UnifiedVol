/**
* CalibratorUtils.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <array>
#include <cstddef>
#include <string_view>
#include <span>

namespace uv
{   
    // Clamp initial guess within upper and lower bounds
    template <std::size_t N>
    void clamp(std::array<double, N>& initGuess,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Warn if upper or lower bounds are touched
    template <std::size_t N>
    void warnBoundsHit(std::span<double> x,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept;

    // Log calibration results 
    template <std::size_t N>
    void logResults(std::span<double> x,
        const std::array<std::string_view, N>& paramNames,
        double sse,
        unsigned iterCount,
        double elapsedMs,
        bool isSuccess) noexcept;
}

#include "CalibratorUtils.inl"