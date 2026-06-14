// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "IO/JSON/Read.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

namespace uv::tests::performance
{
struct Budget
{
    double maxMs{};
};

inline constexpr std::string_view ExamplePipelineBudgetKey{"examplePipeline"};
inline constexpr std::string_view HestonMediumSurfaceBudgetKey{"hestonMediumSurface"};
inline constexpr std::string_view SVISyntheticCalibrationBudgetKey{
    "sviSyntheticCalibration"
};
inline constexpr std::string_view BSplineLargeEvaluationBudgetKey{"bsplineLargeEvaluation"
};

inline constexpr std::array<std::string_view, 4> expectedBudgetKeys()
{
    return {
        ExamplePipelineBudgetKey,
        HestonMediumSurfaceBudgetKey,
        SVISyntheticCalibrationBudgetKey,
        BSplineLargeEvaluationBudgetKey
    };
}

inline Budget readBudget(const std::filesystem::path& path, const std::string_view key)
{
    const auto root = uv::io::json::read(path);
    const std::string keyString{key};
    const double maxMs{root.at("budgets").at(keyString).at("maxMs").asNumber()};
    if (!std::isfinite(maxMs) || maxMs <= 0.0)
        throw std::runtime_error(
            "Performance budget must be positive and finite: " + keyString
        );
    return {.maxMs = maxMs};
}
} // namespace uv::tests::performance
