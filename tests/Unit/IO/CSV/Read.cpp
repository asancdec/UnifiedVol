// SPDX-License-Identifier: Apache-2.0

#include "IO/CSV/Read.hpp"
#include "Base/Errors/Errors.hpp"

#include <gtest/gtest.h>
#include <sstream>

TEST(UnitIOCSVRead, TrimsAndSplitsCommaSeparatedCells)
{
    EXPECT_EQ(uv::io::csv::trimView(" \t 12.5 \n "), "12.5");

    const auto cells = uv::io::csv::splitComma("a,b,");

    ASSERT_EQ(cells.size(), 3U);
    EXPECT_EQ(cells[0], "a");
    EXPECT_EQ(cells[1], "b");
    EXPECT_EQ(cells[2], "");
}

TEST(UnitIOCSVRead, ParsesNumbersAndPercentCells)
{
    EXPECT_DOUBLE_EQ(
        uv::io::csv::parseNumberCellOrThrow<double>(" 12.5% ", "cell", 1, 1),
        0.125
    );
    EXPECT_DOUBLE_EQ(
        uv::io::csv::parseNumberCellOrThrow<double>("-3.25", "cell", 1, 2),
        -3.25
    );
}

TEST(UnitIOCSVRead, ReadsLabeledDenseMatrix)
{
    std::istringstream csv{"maturity,90,100\n"
                           "0.5,20%,21%\n"
                           "1.0,22%,23%\n"};

    const auto dense = uv::io::csv::readLabeledDenseOrThrow<double>(csv, "memory.csv");

    EXPECT_EQ(dense.rows, 2U);
    EXPECT_EQ(dense.cols, 2U);
    EXPECT_DOUBLE_EQ(dense.rowLabels[0], 0.5);
    EXPECT_DOUBLE_EQ(dense.colLabels[1], 100.0);
    EXPECT_DOUBLE_EQ(dense.values[0], 0.20);
    EXPECT_DOUBLE_EQ(dense.values[3], 0.23);
}

TEST(UnitIOCSVRead, RejectsMalformedNumericCells)
{
    EXPECT_THROW(
        uv::io::csv::parseNumberCellOrThrow<double>("abc", "cell", 2, 3),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOCSVRead, RejectsEmptyAndLonelyPercentCells)
{
    EXPECT_THROW(
        uv::io::csv::parseNumberCellOrThrow<double>("   ", "cell", 2, 3),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        uv::io::csv::parseNumberCellOrThrow<double>(" % ", "cell", 2, 3),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOCSVRead, HandlesBlankLinesWhenConfigured)
{
    std::istringstream csv{"maturity,90\n"
                           "\n"
                           "0.5,20%\n"};

    const auto dense = uv::io::csv::readLabeledDenseOrThrow<double>(csv, "memory.csv");

    EXPECT_EQ(dense.rows, 1U);
    EXPECT_EQ(dense.cols, 1U);
    EXPECT_DOUBLE_EQ(dense.values[0], 0.20);
}

TEST(UnitIOCSVRead, RejectsMalformedDenseCsv)
{
    {
        std::istringstream csv{};
        EXPECT_THROW(
            uv::io::csv::readLabeledDenseOrThrow<double>(csv, "empty.csv"),
            uv::errors::UnifiedVolError
        );
    }
    {
        std::istringstream csv{"maturity\n0.5\n"};
        EXPECT_THROW(
            uv::io::csv::readLabeledDenseOrThrow<double>(csv, "bad-header.csv"),
            uv::errors::UnifiedVolError
        );
    }
    {
        std::istringstream csv{"maturity,90,100\n0.5,20%\n"};
        EXPECT_THROW(
            uv::io::csv::readLabeledDenseOrThrow<double>(csv, "short-row.csv"),
            uv::errors::UnifiedVolError
        );
    }
    {
        std::istringstream csv{"maturity,90\n0.5,20%,21%\n"};
        EXPECT_THROW(
            uv::io::csv::readLabeledDenseOrThrow<double>(
                csv,
                "extra-cols.csv",
                {.allowExtraCols = false}
            ),
            uv::errors::UnifiedVolError
        );
    }
    {
        std::istringstream csv{"maturity,90\n"};
        EXPECT_THROW(
            uv::io::csv::readLabeledDenseOrThrow<double>(csv, "no-rows.csv"),
            uv::errors::UnifiedVolError
        );
    }
}
