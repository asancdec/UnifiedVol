// SPDX-License-Identifier: Apache-2.0

#include "Base/Utils/Detail/Log.hpp"
#include "Base/Macros/Require.hpp"

#include <array>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>
#include <string_view>
#include <system_error>

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

    REQUIRE_DIR_CREATED(!ec, logDir);

    std::filesystem::path fullPath = logDir / filename;

    file_.open(fullPath, std::ios::out | std::ios::app);

    fileEnabled_ = file_.is_open();

    REQUIRE_FILE_OPENED(fileEnabled_, fullPath);
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

    std::array<char, 20> timestamp{};
    std::strftime(timestamp.data(), timestamp.size(), "%Y-%m-%d %H:%M:%S", &tm);
    const std::string line{
        std::format("[{}.{:03}][{}] {}\n", timestamp.data(), ms.count(), lvlStr, msg)
    };

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
