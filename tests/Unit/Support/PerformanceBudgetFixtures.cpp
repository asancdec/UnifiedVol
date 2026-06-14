// SPDX-License-Identifier: Apache-2.0

#include "Support/Performance/Budgets.hpp"
#include "Support/TempFile.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(PerformanceBudgetFixtures, CommittedBudgetsHaveValidSchemas)
{
    for (const auto key : uv::tests::performance::expectedBudgetKeys())
    {
        EXPECT_NO_THROW(static_cast<void>(uv::tests::performance::readBudget(
            "tests/Golden/performance_budgets.json",
            key
        ))) << key;
    }
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

TEST(PerformanceBudgetFixtures, RejectsMalformedBudgetShape)
{
    {
        const auto path = uv::tests::writeTempFile(
            "unifiedvol_missing_budgets_object.json",
            R"({"metadata": {}})"
        );
        EXPECT_THROW(
            static_cast<void>(uv::tests::performance::readBudget(
                path,
                uv::tests::performance::ExamplePipelineBudgetKey
            )),
            std::runtime_error
        );
    }
    {
        const auto path = uv::tests::writeTempFile(
            "unifiedvol_missing_max_ms_budget.json",
            R"({"budgets": {"examplePipeline": {}}})"
        );
        EXPECT_THROW(
            static_cast<void>(uv::tests::performance::readBudget(
                path,
                uv::tests::performance::ExamplePipelineBudgetKey
            )),
            std::runtime_error
        );
    }
    {
        const auto path = uv::tests::writeTempFile(
            "unifiedvol_non_numeric_budget.json",
            R"({"budgets": {"examplePipeline": {"maxMs": "fast"}}})"
        );
        EXPECT_THROW(
            static_cast<void>(uv::tests::performance::readBudget(
                path,
                uv::tests::performance::ExamplePipelineBudgetKey
            )),
            std::runtime_error
        );
    }
}

TEST(PerformanceBudgetFixtures, RejectsNonPositiveBudget)
{
    for (const auto maxMs : {0.0, -1.0})
    {
        const auto path = uv::tests::writeTempFile(
            "unifiedvol_bad_performance_budget.json",
            R"({
              "budgets": {
                "examplePipeline": {
                  "maxMs": )" +
                std::to_string(maxMs) + R"(
                }
              }
            })"
        );

        EXPECT_THROW(
            static_cast<void>(uv::tests::performance::readBudget(
                path,
                uv::tests::performance::ExamplePipelineBudgetKey
            )),
            std::runtime_error
        );
    }
}
