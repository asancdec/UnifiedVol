// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <vector>

namespace uv::tests::performance
{
template <typename F>
double bestElapsedMs(F&& f, const std::size_t samples = 3, const std::size_t warmups = 1)
{
    for (std::size_t i = 0; i < warmups; ++i)
        f();

    std::vector<double> timings;
    timings.reserve(samples);
    for (std::size_t i = 0; i < samples; ++i)
    {
        const auto start = std::chrono::steady_clock::now();
        f();
        const auto stop = std::chrono::steady_clock::now();
        timings.emplace_back(
            std::chrono::duration<double, std::milli>{stop - start}.count()
        );
    }

    return *std::min_element(timings.begin(), timings.end());
}
} // namespace uv::tests::performance
