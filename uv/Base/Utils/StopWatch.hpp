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

#include <chrono>
#include <ratio>

#include <string_view>

namespace uv::utils
{

class StopWatch
{
  private:
    std::chrono::high_resolution_clock::time_point startTime_;
    std::chrono::high_resolution_clock::time_point endTime_;
    bool running_{false};

    StopWatch(const StopWatch&) = delete;

    StopWatch& operator=(const StopWatch&) = delete;

  public:
    StopWatch();

    void StartStopWatch() noexcept;

    void StopStopWatch() noexcept;

    void Reset() noexcept;

    template <typename Period = std::milli> double GetTime() const noexcept;

    template <typename Period = std::milli> void LogTime() const noexcept;

    template <typename Period = std::milli>
    void LogTime(std::string_view message) const noexcept;
};
} // namespace uv::utils

#include "Base/Utils/Detail/StopWatch.inl"