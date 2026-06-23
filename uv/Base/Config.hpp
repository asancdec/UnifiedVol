// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

namespace uv
{
struct Config
{
    bool logToConsole{true};
    bool logToFile{true};
    std::string logFile{"calibration.log"};
};

void initialize(const Config& cfg = {});

} // namespace uv
