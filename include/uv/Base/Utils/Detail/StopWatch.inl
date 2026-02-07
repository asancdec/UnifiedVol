// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

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

    UV_INFO(std::format("Elapsed time: {:.6f} {}", dt, unit));
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
        UV_INFO(std::format("{} clocked at: {:.6f} {}", message, dt, unit));
    }
    else
    {
        UV_INFO(std::format("Clocked at: {:.6f} {}", dt, unit));
    }
}
} // namespace uv::utils
