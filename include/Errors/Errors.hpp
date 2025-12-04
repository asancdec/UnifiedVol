/**
* Errors.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <stdexcept>
#include <string>
#include <string_view>
#include <source_location>

namespace uv
{
    //--------------------------------------------------------------------------
    // 1) Domain error code set
    //--------------------------------------------------------------------------
    enum class ErrorCode
    {
        InvalidArgument,
        OutOfRange,
        FileIO,
        NotImplemented,
        DataFormat,
        InvalidState,
        CalibrationError,
        LinearAlgebra,
        Unknown,
    };

    // Converts an ErrorCode to a short string label.
    std::string_view to_string(ErrorCode c) noexcept;

    //--------------------------------------------------------------------------
    // 2) Single exception type carrying code + source location
    //--------------------------------------------------------------------------
    class UnifiedVolError final : public std::runtime_error
    {
    public:
        UnifiedVolError(
            ErrorCode code,
            std::string_view message,
            std::source_location loc = std::source_location::current());

        ErrorCode           code()     const noexcept;
        const std::string&  file()     const noexcept;
        const std::string&  function() const noexcept;
        unsigned            line()     const noexcept;

    private:
        static std::string make_what(ErrorCode code,
            std::string_view msg,
            std::source_location loc);

        ErrorCode   code_;
        std::string file_;
        std::string function_;
        unsigned    line_;
    };

    //--------------------------------------------------------------------------
    // 3) Convenience helpers
    //--------------------------------------------------------------------------
    [[noreturn]] void raise(
        ErrorCode code,
        std::string_view message,
        std::source_location loc = std::source_location::current());

  
 // Guard macro for quick checks (throws UnifiedVolError when cond is false)
#define UV_REQUIRE(cond, code, message) \
        do { if (!(cond)) raise((code), (message)); } while (0)
} // namespace uv

