// SPDX-License-Identifier: Apache-2.0

#include "Base/Utils/StopWatch.hpp"

namespace uv::utils
{
StopWatch::StopWatch() = default;

void StopWatch::StartStopWatch() noexcept
{
    if (!running_)
    {

        running_ = true;

        startTime_ = std::chrono::high_resolution_clock::now();
    }
}

void StopWatch::StopStopWatch() noexcept
{
    if (running_)
    {

        endTime_ = std::chrono::high_resolution_clock::now();

        running_ = false;
    }
}

void StopWatch::Reset() noexcept
{

    running_ = false;
    startTime_ = std::chrono::high_resolution_clock::time_point();
    endTime_ = std::chrono::high_resolution_clock::time_point();
}
} // namespace uv::utils
