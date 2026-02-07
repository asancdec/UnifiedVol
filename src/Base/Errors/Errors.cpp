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

#include "Base/Errors/Errors.hpp"

#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace uv::errors
{

std::string_view to_string(ErrorCode c) noexcept
{
    switch (c)
    {
    case ErrorCode::InvalidArgument:
        return "InvalidArgument";
    case ErrorCode::OutOfRange:
        return "OutOfRange";
    case ErrorCode::FileIO:
        return "FileIO";
    case ErrorCode::NotImplemented:
        return "NotImplemented";
    case ErrorCode::DataFormat:
        return "DataFormat";
    case ErrorCode::CalibrationError:
        return "CalibrationError";
    default:
        return "Unknown";
    }
}

std::string
UnifiedVolError::make_what(ErrorCode code, std::string_view msg, std::source_location loc)
{
    std::ostringstream os;
    os << '[' << to_string(code) << "] " << msg << " @ " << loc.file_name() << ':'
       << loc.line();
    return os.str();
}

UnifiedVolError::UnifiedVolError(
    ErrorCode code,
    std::string_view message,
    std::source_location loc
)
    : std::runtime_error(make_what(code, message, loc)),
      code_(code),
      file_(loc.file_name()),
      function_(loc.function_name()),
      line_(loc.line())
{
}

ErrorCode UnifiedVolError::code() const noexcept
{
    return code_;
}
const std::string& UnifiedVolError::file() const noexcept
{
    return file_;
}
const std::string& UnifiedVolError::function() const noexcept
{
    return function_;
}
unsigned UnifiedVolError::line() const noexcept
{
    return line_;
}

[[noreturn]] void
raise(ErrorCode code, std::string_view message, std::source_location loc)
{
    throw UnifiedVolError(code, message, loc);
}

} // namespace uv::errors