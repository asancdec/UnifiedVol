// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"

#include <source_location>
#include <sstream>
#include <string>
#include <string_view>

namespace uv::errors
{

std::string_view to_string(ErrorCode c) noexcept
{
    using enum ErrorCode;

    switch (c)
    {
    case InvalidArgument:
        return "InvalidArgument";
    case OutOfRange:
        return "OutOfRange";
    case FileIO:
        return "FileIO";
    case NotImplemented:
        return "NotImplemented";
    case DataFormat:
        return "DataFormat";
    case InvalidState:
        return "InvalidState";
    case CalibrationError:
        return "CalibrationError";
    case LinearAlgebra:
        return "LinearAlgebra";
    case Unreachable:
        return "Unreachable";
    case Unknown:
        break;
    }

    return "Unknown";
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
