// SPDX-License-Identifier: Apache-2.0
/*
 * File:        StopWatch.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

#include "Utils/IO/Log.hpp"

#include <format>
#include <type_traits>

namespace uv::utils
{
    template <typename Period>
    Real StopWatch::GetTime() const noexcept
    {
        // If running
        if (running_)
        {
            // If stopwatch is running, calculate the time since start
            auto currentTime_ = std::chrono::high_resolution_clock::now();

            // Return elapsed time
            return std::chrono::duration<Real, Period>(currentTime_ - startTime_).count();
        }
        else
        {
            // If stopwatch is stopped calculate the elapsed time
            return std::chrono::duration<Real, Period>(endTime_ - startTime_).count();
        }
    }

    template <typename Period>
    void StopWatch::LogTime() const noexcept
    {
        const Real dt = GetTime<Period>();

        // Compile-time unit label
        constexpr const char* unit =
            std::is_same_v<Period, std::milli> ? "ms" :
            std::is_same_v<Period, std::micro> ? "us" :
            std::is_same_v<Period, std::nano> ? "ns" :
            "s";

        UV_INFO(std::format("Elapsed time: {:.6f} {}", dt, unit));
    }
}
