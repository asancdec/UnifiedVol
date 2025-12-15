// SPDX-License-Identifier: Apache-2.0
/*
 * File:        StopWatch.cpp
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


#include "Utils/Aux/StopWatch.hpp"
#include <chrono>

namespace uv::utils
{
    StopWatch::StopWatch()   = default;

    void StopWatch::StartStopWatch() noexcept
    {
        if (!running_)
        {
            // Run
            running_ = true;

            // Set startTime to now
            startTime_ = std::chrono::high_resolution_clock::now();
        }
    }

    void StopWatch::StopStopWatch() noexcept
    {
        if (running_)
        {
            // Set the endTime to now
            endTime_ = std::chrono::high_resolution_clock::now();

            // Stop
            running_ = false;
        }
    }

    void StopWatch::Reset() noexcept
    {
        // Reset all private variables
        running_ = false;
        startTime_ = std::chrono::high_resolution_clock::time_point();
        endTime_ = std::chrono::high_resolution_clock::time_point();
    }
}
