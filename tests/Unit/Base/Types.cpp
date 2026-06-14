// SPDX-License-Identifier: Apache-2.0

#include "Base/Types.hpp"

#include <gtest/gtest.h>
#include <span>
#include <vector>

TEST(UnitBaseTypes, ConvertsVectorValues)
{
    const uv::Vector<int> ints{1, 2, 3};

    const auto doubles = uv::convertVector<double>(ints);

    EXPECT_EQ(doubles.size(), ints.size());
    EXPECT_DOUBLE_EQ(doubles[0], 1.0);
    EXPECT_DOUBLE_EQ(doubles[2], 3.0);
}

TEST(UnitBaseTypes, ConvertsSpanValues)
{
    const std::vector<double> doubles{1.25, 2.5};

    const auto floats = uv::convertVector<float>(std::span<const double>{doubles});

    ASSERT_EQ(floats.size(), doubles.size());
    EXPECT_FLOAT_EQ(floats[0], 1.25F);
    EXPECT_FLOAT_EQ(floats[1], 2.5F);
}
