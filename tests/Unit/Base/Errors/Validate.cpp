// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Validate.hpp"
#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>
#include <limits>
#include <optional>
#include <span>
#include <vector>

TEST(UnitBaseErrorsValidate, AcceptsValidScalarAndRangeChecks)
{
    const std::vector<double> increasing{0.0, 1.0, 2.0};
    const std::vector<double> nonDecreasing{0.0, 1.0, 1.0};

    EXPECT_NO_THROW(uv::errors::validate::finite(1.0, "x"));
    EXPECT_NO_THROW(uv::errors::validate::nonNegative(0.0, "x"));
    EXPECT_NO_THROW(uv::errors::validate::positive(1.0, "x"));
    EXPECT_NO_THROW(uv::errors::validate::strictlyIncreasing(increasing, "x"));
    EXPECT_NO_THROW(uv::errors::validate::nonDecreasing(nonDecreasing, "x"));
    EXPECT_NO_THROW(uv::errors::validate::close(1.0, 1.0 + 1e-12, 1e-10, "x"));
}

TEST(UnitBaseErrorsValidate, RejectsInvalidScalarAndRangeChecks)
{
    const std::vector<double> notIncreasing{0.0, 1.0, 1.0};
    const std::vector<double> negative{0.0, -1.0};

    EXPECT_THROW(
        uv::errors::validate::finite(std::numeric_limits<double>::infinity(), "x"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::errors::validate::nonNegative(negative, "x"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(uv::errors::validate::positive(0.0, "x"), uv::errors::UnifiedVolError);
    EXPECT_THROW(
        uv::errors::validate::strictlyIncreasing(notIncreasing, "x"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::errors::validate::close(1.0, 1.1, 1e-3, "x"),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitBaseErrorsValidate, ChecksSizeNullAndOptionalContracts)
{
    const std::vector<double> xs{1.0, 2.0};
    const std::vector<double> ys{3.0, 4.0};
    const std::vector<double> zs{5.0};
    const int value = 3;
    const std::optional<int> setValue{value};
    const std::optional<int> emptyValue{};

    EXPECT_NO_THROW(uv::errors::validate::nonEmpty(xs, "xs"));
    EXPECT_NO_THROW(uv::errors::validate::sameSize(xs, ys, "xs vs ys"));
    EXPECT_NO_THROW(uv::errors::validate::minSize(xs, 2U, "xs"));
    EXPECT_NO_THROW(uv::errors::validate::nonNull(&value, "ptr"));
    EXPECT_NO_THROW(uv::errors::validate::set(setValue, "opt"));

    EXPECT_THROW(
        uv::errors::validate::sameSize(xs, zs, "xs vs zs"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::errors::validate::nonNull(static_cast<const int*>(nullptr), "ptr"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::errors::validate::set(emptyValue, "opt"),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitBaseErrorsValidate, ChecksOrderAgainstThresholds)
{
    const std::vector<double> xs{1.0, 2.0};
    const std::vector<double> lower{0.5, 1.5};
    const std::vector<double> upper{1.5, 2.5};

    EXPECT_NO_THROW(uv::errors::validate::equalOrGreater(
        std::span<const double>{xs},
        std::span<const double>{lower},
        "xs"
    ));
    EXPECT_NO_THROW(uv::errors::validate::equalOrLess(
        std::span<const double>{xs},
        std::span<const double>{upper},
        "xs"
    ));
    EXPECT_THROW(
        uv::errors::validate::greater(1.0, 2.0, "x"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(uv::errors::validate::less(2.0, 1.0, "x"), uv::errors::UnifiedVolError);
}

TEST(UnitBaseErrorsValidate, ChecksDecreasingAndMonotonicContracts)
{
    const std::vector<double> decreasing{3.0, 2.0, 1.0};
    const std::vector<double> increasing{1.0, 2.0, 3.0};
    const std::vector<double> equalPair{1.0, 1.0};
    const std::vector<double> notMonotonic{1.0, 3.0, 2.0};

    EXPECT_NO_THROW(uv::errors::validate::strictlyDecreasing(decreasing, "x"));
    EXPECT_NO_THROW(uv::errors::validate::strictlyMonotonic(decreasing, "x"));
    EXPECT_NO_THROW(uv::errors::validate::strictlyMonotonic(increasing, "x"));
    EXPECT_THROW(
        uv::errors::validate::strictlyMonotonic(equalPair, "x"),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::errors::validate::strictlyMonotonic(notMonotonic, "x"),
        uv::errors::UnifiedVolError
    );
}
