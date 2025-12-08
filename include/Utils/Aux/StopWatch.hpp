// SPDX-License-Identifier: Apache-2.0
/*
 * File:        StopWatch.hpp
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

#include "Utils/Types.hpp"

#include <chrono>

namespace uv::utils
{
    class StopWatch
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------
        std::chrono::high_resolution_clock::time_point startTime_;
        std::chrono::high_resolution_clock::time_point endTime_;
        bool running_;

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------   
        // Private copy constructor to prevent copying
        StopWatch(const StopWatch&);

        // Private copy assignment operator to prevent assignment
        StopWatch& operator=(const StopWatch&);

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------   
        StopWatch();

        //--------------------------------------------------------------------------
        // Time
        //-------------------------------------------------------------------------   

        // Starts the stopwatch if it is not already running
        void StartStopWatch() noexcept;

        // Stops the stopwatch if it is running
        void StopStopWatch() noexcept;

        // Resets the stopwatch clearing the running state and time points
        void Reset() noexcept;

        // Returns the elapsed time
        template <typename Period = std::ratio<1>>
        Real GetTime() const noexcept;

        // Log time
        template <typename Period = std::ratio<1>>
        void LogTime() const noexcept;
    };
}

#include "StopWatch.inl"

