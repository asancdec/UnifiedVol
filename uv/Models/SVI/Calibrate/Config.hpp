// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <string_view>

namespace uv::models::svi
{

struct Config
{
    double objectiveTol{1e-12};
    unsigned int maxEval{10000};
    bool verbose{true};
    bool printParams{false};
};

namespace detail
{
inline constexpr std::array<std::string_view, 4> paramNames{"b", "rho", "m", "sigma"};
}

} // namespace uv::models::svi