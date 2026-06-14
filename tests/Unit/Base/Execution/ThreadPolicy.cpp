// SPDX-License-Identifier: Apache-2.0

#include "Base/Execution/ThreadPolicy.hpp"
#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>

TEST(UnitBaseExecutionThreadPolicy, PositiveRequestsReturnAtLeastOneThread)
{
    EXPECT_GE(uv::execution::requestThreads(1), 1);
}

TEST(UnitBaseExecutionThreadPolicy, NegativeRequestsCountBackFromAvailableThreads)
{
    EXPECT_GE(uv::execution::requestThreads(-1), 1);
}

TEST(UnitBaseExecutionThreadPolicy, RejectsZeroThreads)
{
    EXPECT_THROW(uv::execution::requestThreads(0), uv::errors::UnifiedVolError);
}
