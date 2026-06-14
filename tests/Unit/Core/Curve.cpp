// SPDX-License-Identifier: Apache-2.0

#include "Core/Curve.hpp"
#include "Base/Errors/Errors.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(CoreCurve, ReturnsExactDiscountFactorsAtKnownMaturities)
{
    const std::vector<double> maturities{0.5, 1.0, 2.0};
    const uv::core::Curve<double> curve{0.03, maturities};

    EXPECT_NEAR(curve.interpolateDF(0.5), std::exp(-0.03 * 0.5), 1e-15);
    EXPECT_NEAR(curve.interpolateDF(1.0), std::exp(-0.03), 1e-15);
    EXPECT_NEAR(curve.interpolateDF(2.0), std::exp(-0.06), 1e-15);
}

TEST(CoreCurve, VectorizedInterpolationMatchesScalarInterpolation)
{
    const std::vector<double> maturities{0.25, 1.0, 3.0};
    const uv::core::Curve<double> curve{0.04, maturities};

    const auto dfs = curve.interpolateDF(maturities);

    ASSERT_EQ(dfs.size(), maturities.size());
    for (std::size_t i = 0; i < maturities.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(dfs[i], curve.interpolateDF(maturities[i]));
    }
}

TEST(CoreCurve, RejectsInvalidInputGrid)
{
    const std::vector<double> unsorted{1.0, 0.5};

    EXPECT_THROW((uv::core::Curve<double>{0.03, unsorted}), uv::errors::UnifiedVolError);
}
