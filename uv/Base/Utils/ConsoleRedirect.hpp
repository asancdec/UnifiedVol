// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Macros/Inform.hpp"

#include <iostream>
#include <sstream>

namespace uv::utils
{

struct ConsoleRedirect
{
    explicit ConsoleRedirect()
        : oldBuf_(std::cout.rdbuf(buffer_.rdbuf()))
    {
    }

    ~ConsoleRedirect()
    {
        std::cout.rdbuf(oldBuf_);
        INFO("\n" + buffer_.str());
    }

  private:
    std::stringstream buffer_;
    std::streambuf* oldBuf_{};
};

} // namespace uv::utils
