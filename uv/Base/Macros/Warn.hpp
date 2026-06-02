// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Utils/Detail/Log.hpp"

#define WARN(cond, msg)                                                                  \
    do                                                                                   \
    {                                                                                    \
        if (cond)                                                                        \
        {                                                                                \
            ::uv::utils::Log::instance().log(::uv::utils::Level::Warn, (msg));           \
        }                                                                                \
    } while (0)
