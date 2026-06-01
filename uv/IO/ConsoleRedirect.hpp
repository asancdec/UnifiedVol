// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Macros/Inform.hpp"

#include <iostream>
#include <sstream>

namespace uv::io
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
        UV_INFO("\n" + buffer_.str());
    }

  private:
    std::stringstream buffer_;
    std::streambuf* oldBuf_{};
};
} // namespace uv::io