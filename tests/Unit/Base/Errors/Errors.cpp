// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>

namespace uv::tests::unit::base_errors::detail
{
enum class DummyEnum
{
    Known = 1,
    Unknown = 7
};
} // namespace uv::tests::unit::base_errors::detail

using namespace uv::tests::unit::base_errors::detail;

TEST(UnitBaseErrors, UnreachableEnumRaisesUnifiedVolError)
{
    EXPECT_THROW(
        uv::errors::unreachableEnum(DummyEnum::Unknown, "DummyEnum"),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitBaseErrors, ErrorCodesHaveStableNames)
{
    using enum uv::errors::ErrorCode;

    EXPECT_EQ(uv::errors::to_string(InvalidState), "InvalidState");
    EXPECT_EQ(uv::errors::to_string(LinearAlgebra), "LinearAlgebra");
    EXPECT_EQ(uv::errors::to_string(Unreachable), "Unreachable");
    EXPECT_EQ(uv::errors::to_string(Unknown), "Unknown");
}
