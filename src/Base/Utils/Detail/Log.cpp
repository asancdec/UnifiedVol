// SPDX-License-Identifier: Apache-2.0
/*
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

#include <Base/Macros/Require.hpp>
#include <Base/Utils/Detail/Log.hpp>

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

namespace uv::utils
{

using errors::ErrorCode;

Log& Log::instance()
{
    static Log g;
    return g;
}

Log::Log() = default;

void Log::setFile(std::string_view filename)
{
    std::filesystem::path root = std::filesystem::current_path();

    while (!std::filesystem::exists(root / "CMakeLists.txt") && root.has_parent_path())
        root = root.parent_path();

    std::filesystem::path logDir = root / "logs";
    std::error_code ec;

    std::filesystem::create_directories(logDir, ec);

    UV_REQUIRE_DIR_CREATED(!ec, logDir);

    std::filesystem::path fullPath = logDir / filename;

    file_.open(fullPath, std::ios::out | std::ios::app);

    fileEnabled_ = file_.is_open();

    UV_REQUIRE_FILE_OPENED(fileEnabled_, fullPath);
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

    if (consoleEnabled_)
    {
        std::cout << line << std::flush;
    }

    if (fileEnabled_)
    {
        file_ << line;
        file_.flush();
    }
}

} // namespace uv::utils