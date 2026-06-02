// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Utils/Detail/Log.hpp"

#define INFO(msg)                                                                        \
    do                                                                                   \
    {                                                                                    \
        ::uv::utils::Log::instance().log(::uv::utils::Level::Info, (msg));               \
    } while (0)
