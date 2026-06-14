// SPDX-License-Identifier: Apache-2.0

#include "Budgets.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace
{
std::filesystem::path
writeBudgetFixture(const std::string& name, const std::string& contents)
{
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out{path};
    out << contents;
    return path;
}
} // namespace

TEST(PerformanceBudgetFixtures, CommittedBudgetsHaveValidSchemas)
{
    EXPECT_NO_THROW(static_cast<void>(uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        "examplePipeline"
    )));
    EXPECT_NO_THROW(static_cast<void>(uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        "hestonMediumSurface"
    )));
    EXPECT_NO_THROW(static_cast<void>(uv::tests::performance::readBudget(
        "tests/Golden/performance_budgets.json",
        "sviSyntheticCalibration"
    )));
}

TEST(PerformanceBudgetFixtures, RejectsMissingBudgetKey)
{
    EXPECT_THROW(
        static_cast<void>(uv::tests::performance::readBudget(
            "tests/Golden/performance_budgets.json",
            "missingBudget"
        )),
        std::runtime_error
    );
}

TEST(PerformanceBudgetFixtures, RejectsNonPositiveBudget)
{
    const auto path = writeBudgetFixture(
        "unifiedvol_bad_performance_budget.json",
        R"({
          "budgets": {
            "examplePipeline": {
              "maxMs": 0.0
            }
          }
        })"
    );

    EXPECT_THROW(
        static_cast<void>(uv::tests::performance::readBudget(path, "examplePipeline")),
        std::runtime_error
    );
}
