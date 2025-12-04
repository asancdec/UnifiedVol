/**
* StopWatch.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <chrono>

namespace uv
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
        double GetTime() const noexcept;

        // Log time
        template <typename Period = std::ratio<1>>
        void LogTime() const noexcept;
    };
}

#include "StopWatch.inl"

