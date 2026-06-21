// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace uv::opt::nlopt
{

template <std::size_t N> struct Config
{
    double eps{1e-12};
    double tol{};
    double ftolRel{};
    unsigned int maxEval{};
    bool verbose{};
    std::array<std::string_view, N> paramNames;
};
} // namespace uv::opt::nlopt
