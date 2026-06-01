// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Errors/Errors.hpp"

#include <source_location>

#define UV_UNREACHABLE(EnumType, value)                                                  \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::unreachableEnum(                                                   \
            (value),                                                                     \
            #EnumType,                                                                   \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)
