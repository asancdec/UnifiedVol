// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Errors/Errors.hpp"

#include <source_location>

#define NOT_IMPLEMENTED(msg)                                                             \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::raise(                                                             \
            ::uv::errors::ErrorCode::NotImplemented,                                     \
            (msg),                                                                       \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)
