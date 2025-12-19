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

#include <chrono>

namespace uv::utils
{
    /**
     * @brief Lightweight RAII-style stopwatch for wall-clock timing.
     *
     * Provides simple start/stop/reset semantics and supports querying
     * elapsed time in arbitrary chrono periods.
     *
     * Typical use cases:
     * - timing calibration routines
     * - profiling numerical kernels
     * - logging elapsed wall-clock time
     *
     * @note This class measures wall-clock time, not CPU time.
     * @note Not thread-safe.
     */

    class StopWatch
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        std::chrono::high_resolution_clock::time_point startTime_; // Start timestamp
        std::chrono::high_resolution_clock::time_point endTime_;   // End timestamp
        bool running_{false};                                             // Whether the stopwatch is running

        //--------------------------------------------------------------------------
        // Copy prevention
        //--------------------------------------------------------------------------

        /// Deleted copy constructor (StopWatch is non-copyable).
        StopWatch(const StopWatch&) = delete;

        /// Deleted copy assignment operator (StopWatch is non-copyable).
        StopWatch& operator=(const StopWatch&) = delete;

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------

        /**
         * @brief Construct a stopped stopwatch with zero elapsed time.
         */
        StopWatch();

        //--------------------------------------------------------------------------
        // Time control
        //--------------------------------------------------------------------------

        /**
         * @brief Start the stopwatch.
         *
         * Has no effect if the stopwatch is already running.
         */
        void StartStopWatch() noexcept;

        /**
         * @brief Stop the stopwatch.
         *
         * Has no effect if the stopwatch is not running.
         */
        void StopStopWatch() noexcept;

        /**
         * @brief Reset the stopwatch.
         *
         * Clears the running state and stored timestamps.
         */
        void Reset() noexcept;

        /**
         * @brief Get the elapsed time.
         *
         * If the stopwatch is currently running, returns the time elapsed
         * since the last call to @ref StartStopWatch.
         * Otherwise, returns the time between the last start/stop pair.
         *
         * @tparam Period Chrono period for the returned duration
         *                (e.g. std::milli, std::micro).
         * @return Elapsed time expressed in the given period.
         */
        template <typename Period = std::ratio<1>>
        double GetTime() const noexcept;

        /**
         * @brief Log the elapsed time using the UnifiedVol logger.
         *
         * Formats the elapsed time with a unit label derived from @tparam Period
         * and emits an informational log message.
         *
         * @tparam Period Chrono period for the logged duration.
         */
        template <typename Period = std::ratio<1>>
        void LogTime() const noexcept;
    };
} // namespace uv::utils

#include "StopWatch.inl"