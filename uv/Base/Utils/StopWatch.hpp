// SPDX-License-Identifier: Apache-2.0

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