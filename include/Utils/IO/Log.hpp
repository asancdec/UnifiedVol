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
    /**
     * @brief Log severity level.
     */
    enum class Level
    {
        Info,   ///< Informational message
        Warn    ///< Warning message
    };

    /**
     * @brief Central logging facility for UnifiedVol.
     *
     * Implements a lightweight singleton logger supporting:
     * - console output (stdout)
     * - file output under <project_root>/logs/
     *
     * Logging is timestamped with millisecond precision and formatted
     * consistently across the library.
     *
     * The logger is configured once at startup (typically in @c main()).
     *
     * @note This class is not thread-safe.
     * @note File output is optional and disabled by default.
     */
    class Log
    {
    public:
        /**
         * @brief Retrieve the global logger instance.
         *
         * Uses the Meyers singleton pattern.
         *
         * @return Reference to the global @ref Log instance.
         */
        static Log& instance();

        /**
         * @brief Enable file logging.
         *
         * Opens a log file located at:
         *   <project_root>/logs/<path>
         *
         * The project root is detected by walking up the directory tree until
         * a CMakeLists.txt file is found.
         *
         * @param path Log file name (relative, without directories).
         *
         * @throws UnifiedVolError if the directory or file cannot be created.
         */
        void setFile(std::string_view path);

        /**
         * @brief Enable or disable console logging.
         *
         * @param enabled If true, log messages are printed to std::cout.
         */
        void enableConsole(bool enabled) noexcept;

        /**
         * @brief Core logging API.
         *
         * Formats and writes a log entry with timestamp and severity level.
         * The message is written to:
         * - std::cout (if console logging is enabled)
         * - the log file (if file logging is enabled)
         *
         * @param lvl Log severity level.
         * @param msg Message to log.
         */
        void log(Level lvl, std::string_view msg);

    private:
        /**
         * @brief Private constructor to enforce singleton pattern.
         */
        Log();

        std::ofstream file_;     ///< Log file stream
        bool fileEnabled_{ false };///< Whether file logging is enabled
        bool consoleEnabled_{ true }; ///< Whether console logging is enabled
    };

    // -------------------------------------------------------------------------
    // Logging macros
    // -------------------------------------------------------------------------

    /**
     * @brief Enable file logging.
     *
     * Typically called once at program startup.
     *
     * @param pathstr Log file name (relative to <project_root>/logs).
     */
#define UV_LOG_TO_FILE(pathstr) \
        Log::instance().setFile((pathstr))

     /**
      * @brief Enable or disable console logging.
      *
      * @param enabled Boolean flag.
      */
#define UV_LOG_CONSOLE(enabled) \
        Log::instance().enableConsole((enabled))

      /**
       * @brief Emit an informational log message.
       *
       * @param msg Message to log.
       */
#define UV_INFO(msg) \
        Log::instance().log(Level::Info, (msg))

       /**
        * @brief Emit a warning log message conditionally.
        *
        * The message is logged only if the condition evaluates to true.
        *
        * @param cond Condition to test.
        * @param msg  Message to log if condition holds.
        */
#define UV_WARN(cond, msg) \
        do { \
            if (cond) \
                Log::instance().log(Level::Warn, (msg)); \
        } while (0)

} // namespace uv


