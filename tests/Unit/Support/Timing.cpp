// SPDX-License-Identifier: Apache-2.0

#include "Support/Performance/Timing.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(PerformanceTiming, RejectsZeroSamples)
{
    EXPECT_THROW(
        static_cast<void>(uv::tests::performance::bestElapsedMs([] {}, 0, 0)),
        std::invalid_argument
    );
}

TEST(PerformanceTiming, InvokesWarmupsBeforeMeasuredSamples)
{
    int calls{};

    const double elapsedMs = uv::tests::performance::bestElapsedMs(
        [&]
        {
            ++calls;
        },
        3,
        2
    );

    EXPECT_EQ(calls, 5);
    EXPECT_GE(elapsedMs, 0.0);
}
