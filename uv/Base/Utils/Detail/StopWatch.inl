// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Macros/Inform.hpp"

#include <format>
#include <type_traits>

namespace uv::utils
{
template <typename Period> double StopWatch::GetTime() const noexcept
{

    if (running_)
    {

        auto currentTime_ = std::chrono::high_resolution_clock::now();

        return std::chrono::duration<double, Period>(currentTime_ - startTime_).count();
    }
    else
    {

        return std::chrono::duration<double, Period>(endTime_ - startTime_).count();
    }
}

template <typename Period> void StopWatch::LogTime() const noexcept
{
    const double dt = GetTime<Period>();

    constexpr const char* unit = std::is_same_v<Period, std::milli>   ? "ms"
                                 : std::is_same_v<Period, std::micro> ? "us"
                                 : std::is_same_v<Period, std::nano>  ? "ns"
                                                                      : "s";

    INFO(std::format("Elapsed time: {:.6f} {}", dt, unit));
}

template <typename Period>
void StopWatch::LogTime(std::string_view message) const noexcept
{
    const double dt = GetTime<Period>();

    constexpr const char* unit = std::is_same_v<Period, std::milli>   ? "ms"
                                 : std::is_same_v<Period, std::micro> ? "us"
                                 : std::is_same_v<Period, std::nano>  ? "ns"
                                                                      : "s";

    if (!message.empty())
    {
        INFO(std::format("{} clocked at: {:.6f} {}", message, dt, unit));
    }
    else
    {
        INFO(std::format("Clocked at: {:.6f} {}", dt, unit));
    }
}
} // namespace uv::utils
