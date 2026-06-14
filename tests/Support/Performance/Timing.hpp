// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Utils/StopWatch.hpp"

#include <algorithm>
#include <cstddef>
#include <ratio>
#include <stdexcept>
#include <vector>

namespace uv::tests::performance
{
template <typename F>
double bestElapsedMs(F&& f, const std::size_t samples = 3, const std::size_t warmups = 1)
{
    if (samples == 0)
        throw std::invalid_argument("Performance timing requires at least one sample");

    for (std::size_t i = 0; i < warmups; ++i)
        f();

    std::vector<double> timings;
    timings.reserve(samples);
    for (std::size_t i = 0; i < samples; ++i)
    {
        uv::utils::StopWatch watch;
        watch.StartStopWatch();
        f();
        watch.StopStopWatch();
        timings.emplace_back(watch.GetTime<std::milli>());
    }

    return *std::min_element(timings.begin(), timings.end());
}
} // namespace uv::tests::performance
