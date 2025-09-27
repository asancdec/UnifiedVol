// Log.cpp
// Author: Alvaro Sanchez de Carlos
// Date: 09/14/2025

#include "Utils/Log.hpp"

#include <iostream>

namespace uv 
{

    void warn(std::string_view msg, std::source_location loc)
    {
        std::cerr << "[Warning] " << msg
            << " @ " << loc.file_name() << ':' << loc.line() << '\n';
    }

} // namespace uv