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

    void setFile(std::string_view path);

    void enableConsole(bool enabled) noexcept;

    void log(Level lvl, std::string_view msg);

  private:
    Log();

    std::ofstream file_;
    bool fileEnabled_{false};
    bool consoleEnabled_{true};
};
} // namespace uv::utils
