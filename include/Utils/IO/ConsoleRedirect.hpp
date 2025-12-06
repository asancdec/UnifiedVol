/**
* ConsoleRedirect.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/IO/Log.hpp"

#include <iostream>   
#include <sstream>     
#include <string>       

namespace uv::utils
{

    struct ConsoleRedirect
    {
        explicit ConsoleRedirect()
        {
            oldBuf_ = std::cout.rdbuf(buffer_.rdbuf());
        }

        ~ConsoleRedirect()
        {
            std::cout.rdbuf(oldBuf_);
            UV_INFO("\n" +  buffer_.str());
        }

        private:
            std::stringstream buffer_;
            std::streambuf* oldBuf_{};
    };}