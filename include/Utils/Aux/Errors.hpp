// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Errors.hpp
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

#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>

#pragma once

#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>

namespace uv
{
//--------------------------------------------------------------------------
// 1) Domain error code set
//--------------------------------------------------------------------------

/**
 * @brief UnifiedVol domain error codes.
 *
 * Identifies broad classes of failures (validation, I/O, calibration, etc.)
 * and is carried by @ref UnifiedVolError for structured error handling.
 *
 * These codes are intended to be stable and used for:
 * - diagnostics and logging
 * - mapping failures to exit codes in executables
 * - testing expected failure modes
 */
enum class ErrorCode
{
    InvalidArgument,  // Invalid argument passed by caller
    OutOfRange,       // Index/value outside valid range
    FileIO,           // File or filesystem I/O error
    NotImplemented,   // Feature not implemented
    DataFormat,       // Malformed or inconsistent input data
    InvalidState,     // Object/system in invalid state for operation
    CalibrationError, // Calibration/optimisation failed or unstable
    LinearAlgebra,    // Linear algebra failure (singular, factorisation, etc.)
    Unknown,          // Fallback / unknown error
};

/**
 * @brief Convert an ErrorCode to a short string label.
 *
 * @param c Error code.
 * @return String label (e.g. "FileIO", "CalibrationError").
 */
std::string_view to_string(ErrorCode c) noexcept;

//--------------------------------------------------------------------------
// 2) Single exception type carrying code + source location
//--------------------------------------------------------------------------

/**
 * @brief UnifiedVol exception type carrying an error code and source location.
 *
 * Extends std::runtime_error with:
 * - a domain-specific @ref ErrorCode
 * - file/function/line captured via std::source_location
 *
 * @note The @c what() string is preformatted and includes the code + message
 *       and source location (file:line).
 */
class UnifiedVolError final : public std::runtime_error
{
  public:
    /**
     * @brief Construct an exception with code, message, and source location.
     *
     * @param code Domain error code.
     * @param message Human-readable error message.
     * @param loc Source location of the throw site (defaults to current()).
     */
    UnifiedVolError(
        ErrorCode code,
        std::string_view message,
        std::source_location loc = std::source_location::current()
    );

    /** @brief Get the associated domain error code. */
    ErrorCode code() const noexcept;

    /** @brief Source file name where the error originated. */
    const std::string& file() const noexcept;

    /** @brief Function name where the error originated. */
    const std::string& function() const noexcept;

    /** @brief Source line number where the error originated. */
    unsigned line() const noexcept;

  private:
    /**
     * @brief Build the std::runtime_error::what() string.
     *
     * Formats: "[<Code>] <message> @ <file>:<line>"
     */
    static std::string
    make_what(ErrorCode code, std::string_view msg, std::source_location loc);

    ErrorCode code_;       // Domain error code
    std::string file_;     // Source file name
    std::string function_; // Source function name
    unsigned line_;        // Source line number
};

//--------------------------------------------------------------------------
// 3) Convenience helpers
//--------------------------------------------------------------------------

/**
 * @brief Throw a @ref UnifiedVolError.
 *
 * @param code Domain error code.
 * @param message Human-readable error message.
 * @param loc Source location (defaults to current()).
 *
 * @throws UnifiedVolError Always.
 */
[[noreturn]] void raise(
    ErrorCode code,
    std::string_view message,
    std::source_location loc = std::source_location::current()
);

//--------------------------------------------------------------------------
// 4) Guard macro
//--------------------------------------------------------------------------

/**
 * @brief Guard macro for quick checks.
 *
 * Throws a @ref UnifiedVolError with the given @ref ErrorCode and message
 * when the condition evaluates to false.
 *
 * @param cond Condition to validate.
 * @param code Error code to throw on failure.
 * @param message Error message to throw on failure.
 */
#define UV_REQUIRE(cond, code, message)                                                  \
    do                                                                                   \
    {                                                                                    \
        if (!(cond))                                                                     \
            raise((code), (message));                                                    \
    } while (0)

} // namespace uv
