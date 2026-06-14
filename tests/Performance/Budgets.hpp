// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Support/Json.hpp"

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace uv::tests::performance
{
struct Budget
{
    double maxMs{};
};

inline Budget readBudget(const std::filesystem::path& path, const std::string& key)
{
    const auto root = json::read(path);
    const double maxMs{root.at("budgets").at(key).at("maxMs").asNumber()};
    if (!std::isfinite(maxMs) || maxMs <= 0.0)
        throw std::runtime_error(
            "Performance budget must be positive and finite: " + key
        );
    return {.maxMs = maxMs};
}
} // namespace uv::tests::performance
