/**
* StopWatch.inl
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Log.hpp"

#include <format>
#include <type_traits>

namespace uv
{
    template <typename Period>
    double StopWatch::GetTime() const noexcept
    {
        // If running
        if (running_)
        {
            // If stopwatch is running, calculate the time since start
            auto currentTime_ = std::chrono::high_resolution_clock::now();

            // Return elapsed time
            return std::chrono::duration<double, Period>(currentTime_ - startTime_).count();
        }
        else
        {
            // If stopwatch is stopped calculate the elapsed time
            return std::chrono::duration<double, Period>(endTime_ - startTime_).count();
        }
    }

    template <typename Period>
    void StopWatch::LogTime() const noexcept
    {
        const double dt = GetTime<Period>();

        // Compile-time unit label
        constexpr const char* unit =
            std::is_same_v<Period, std::milli> ? "ms" :
            std::is_same_v<Period, std::micro> ? "µs" :
            std::is_same_v<Period, std::nano> ? "ns" :
            "s";

        UV_INFO(std::format("Elapsed time: {:.6f} {}", dt, unit));
    }
}
