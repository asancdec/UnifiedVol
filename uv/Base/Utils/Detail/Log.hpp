// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <fstream>
#include <string_view>

namespace uv::utils
{

enum class Level
{
    Info,
    Warn
};

class Log
{
  public:
    static Log& instance();

    void setFile(std::string_view filename);

    void enableConsole(bool enabled) noexcept;

    void log(Level lvl, std::string_view msg);

  private:
    Log();

    std::ofstream file_;
    bool fileEnabled_{false};
    bool consoleEnabled_{true};
};
} // namespace uv::utils
