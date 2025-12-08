// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Log.hpp
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


#pragma once
#include <fstream>
#include <string_view>

namespace uv
{
    enum class Level { Info, Warn };

    class Log
    {
    public:
        // Singleton instance
        static Log& instance();

        // Enable file output (writes to <project_root>/logs/<path>).
        void setFile(std::string_view path);

        // Enable or disable console output
        void enableConsole(bool enabled) noexcept;

        // Core logging API (prints to console; also to file if enabled).
        void log(Level lvl, std::string_view msg);

    private:
        Log(); // private constructor to enforce singleton pattern
        std::ofstream file_;
        bool fileEnabled_{ false };
        bool consoleEnabled_{ true };
    };

    // ---------------------------------------------------------------------------
    // Macros
    // ---------------------------------------------------------------------------

    // Configure file output once at startup (e.g. in main()).
    #define UV_LOG_TO_FILE(pathstr) \
        Log::instance().setFile((pathstr))

    // Enable or disable console logging
    #define UV_LOG_CONSOLE(enabled) \
        Log::instance().enableConsole((enabled))

    // Info-level message (always printed)
    #define UV_INFO(msg) \
        Log::instance().log(Level::Info, (msg))

    // Warning message, printed only if condition is true
    #define UV_WARN(cond, msg) \
        do { \
            if (cond) \
                Log::instance().log(Level::Warn, (msg)); \
        } while (0)
} // namespace uv


