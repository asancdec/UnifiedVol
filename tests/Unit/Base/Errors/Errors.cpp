// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>

namespace
{
enum class DummyEnum
{
    Known = 1,
    Unknown = 7
};
} // namespace

TEST(UnitBaseErrors, UnreachableEnumRaisesUnifiedVolError)
{
    EXPECT_THROW(
        uv::errors::unreachableEnum(DummyEnum::Unknown, "DummyEnum"),
        uv::errors::UnifiedVolError
    );
}
