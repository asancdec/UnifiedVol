// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Log.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#include "Utils/IO/Log.hpp"
#include "Utils/Aux/Errors.hpp"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <time.h>

namespace uv
{

Log& Log::instance()
{
    static Log g;
    return g;
}

Log::Log() = default;

void Log::setFile(std::string_view filename)
{
    namespace fs = std::filesystem;

    fs::path root = fs::current_path();

    while (!fs::exists(root / "CMakeLists.txt") && root.has_parent_path())
        root = root.parent_path();

    fs::path logDir = root / "logs";
    std::error_code ec;
    fs::create_directories(logDir, ec);
    UV_REQUIRE(
        !ec,
        ErrorCode::FileIO,
        "Failed to create logs directory: " + logDir.string()
    );

    fs::path fullPath = logDir / filename;

    file_.open(fullPath, std::ios::out | std::ios::app);
    fileEnabled_ = file_.is_open();

    UV_REQUIRE(
        fileEnabled_,
        ErrorCode::FileIO,
        "Unable to open log file: " + fullPath.string()
    );
}

void Log::enableConsole(bool enabled) noexcept
{
    consoleEnabled_ = enabled;
}

void Log::log(Level lvl, std::string_view msg)
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const auto tt = system_clock::to_time_t(now);
    const auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &tt);
#else
    localtime_r(&tt, &tm);
#endif

    const char* lvlStr = (lvl == Level::Info) ? "INFO" : "WARN";

    std::ostringstream oss;
    oss << '[' << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << '.' << std::setw(3)
        << std::setfill('0') << ms.count() << ']' << '[' << lvlStr << "] " << msg << '\n';

    const std::string line = oss.str();

    // Console output (toggle)
    if (consoleEnabled_)
    {
        std::cout << line << std::flush;
    }

    // File sink
    if (fileEnabled_)
    {
        file_ << line;
        file_.flush();
    }
}

} // namespace uv