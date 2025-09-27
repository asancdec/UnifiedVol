// Log.hpp
// Author: Alvaro Sanchez de Carlos
// Date: 09/14/2025

#ifndef UV_LOG_HPP
#define UV_LOG_HPP

#include <source_location>
#include <string_view>

namespace uv {

    // Print a warning with file:line (no throw)
    void warn(std::string_view msg,
        std::source_location loc = std::source_location::current());

} // namespace uv

// Warn if NOT(cond). Just prints; never throws.
#define UV_WARN(cond, message) \
    do { if (!(cond)) ::uv::warn((message)); } while (0)

#endif // UV_LOG_HPP