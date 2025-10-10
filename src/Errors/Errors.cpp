/**
* Errors.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Errors/Errors.hpp"

#include <sstream>  

namespace uv
{
    //--------------------------------------------------------------------------
    // ErrorCode → string
    //--------------------------------------------------------------------------
    std::string_view to_string(ErrorCode c) noexcept
    {
        using E = ErrorCode;
        switch (c) {
        case E::InvalidArgument:  return "InvalidArgument";
        case E::OutOfRange:       return "OutOfRange";
        case E::FileIO:           return "FileIO";
        case E::NotImplemented:   return "NotImplemented";
        case E::DataFormat:       return "DataFormat";
        case E::CalibrationError: return "CalibrationError";
        default:                  return "Unknown";
        }
    }

    //--------------------------------------------------------------------------
    // UnifiedVolError internals
    //--------------------------------------------------------------------------
    std::string UnifiedVolError::make_what(
        ErrorCode code,
        std::string_view msg,
        std::source_location loc)
    {
        std::ostringstream os;
        os << '[' << to_string(code) << "] " << msg
            << " @ " << loc.file_name() << ':' << loc.line();
        return os.str();
    }

    UnifiedVolError::UnifiedVolError(
        ErrorCode code,
        std::string_view message,
        std::source_location loc)
        : std::runtime_error(make_what(code, message, loc))
        , code_(code)
        , file_(loc.file_name())
        , function_(loc.function_name())
        , line_(loc.line())
    {
    }

    ErrorCode UnifiedVolError::code() const noexcept { return code_; }
    const std::string& UnifiedVolError::file() const noexcept { return file_; }
    const std::string& UnifiedVolError::function() const noexcept { return function_; }
    unsigned UnifiedVolError::line() const noexcept { return line_; }

    //--------------------------------------------------------------------------
    // raise helper
    //--------------------------------------------------------------------------
    [[noreturn]] void raise(
        ErrorCode code,
        std::string_view message,
        std::source_location loc)
    {
        throw UnifiedVolError(code, message, loc);
    }

} // namespace uv
