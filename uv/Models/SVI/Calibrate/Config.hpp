// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace uv::models::svi
{

struct Config
{
    double objectiveTol{1e-12};
    unsigned int maxEval{10000};
    bool verbose{true};
    bool printParams{false};
};

} // namespace uv::models::svi
