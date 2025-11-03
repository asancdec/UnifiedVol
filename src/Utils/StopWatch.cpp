/**
* StopWatch.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/StopWatch/StopWatch.hpp"

namespace uv
{
    StopWatch::StopWatch() : startTime_(), endTime_(), running_(false) {}

    void StopWatch::StartStopWatch() noexcept
    {
        if (!running_)
        {
            // Run
            running_ = true;

            // Set startTime to now
            startTime_ = ::std::chrono::high_resolution_clock::now();
        }
    }

    void StopWatch::StopStopWatch() noexcept
    {
        if (running_)
        {
            // Set the endTime to now
            endTime_ = ::std::chrono::high_resolution_clock::now();

            // Stop
            running_ = false;
        }
    }

    void StopWatch::Reset() noexcept
    {
        // Reset all private variables
        running_ = false;
        startTime_ = ::std::chrono::high_resolution_clock::time_point();
        endTime_ = ::std::chrono::high_resolution_clock::time_point();
    }
}
