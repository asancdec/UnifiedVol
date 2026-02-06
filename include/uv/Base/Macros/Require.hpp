// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#pragma once

#include <Base/Errors/Errors.hpp>
#include <Base/Errors/Validate.hpp>

#include <source_location>

#define UV_REQUIRE(cond, code, message)                                                  \
    do                                                                                   \
    {                                                                                    \
        if (!cond)                                                                       \
        {                                                                                \
            ::uv::errors::raise((code), (message), std::source_location::current());     \
        }                                                                                \
    } while (0)

#define UV_REQUIRE_FINITE(x)                                                             \
    ::uv::errors::validateFinite((x), #x, std::source_location::current())

#define UV_REQUIRE_NON_NEGATIVE(x)                                                       \
    ::uv::errors::validateNonNegative((x), #x, std::source_location::current())

#define UV_REQUIRE_POSITIVE(x)                                                           \
    ::uv::errors::validatePositive((x), #x, std::source_location::current())

#define UV_REQUIRE_EQUAL_OR_GREATER(x, threshold)                                        \
    ::uv::errors::validateEqualOrGreater(                                                \
        (threshold),                                                                     \
        (x),                                                                             \
        #x,                                                                              \
        std::source_location::current()                                                  \
    )

#define UV_REQUIRE_GREATER(x, threshold)                                                 \
    ::uv::errors::validateGreater((threshold), (x), #x, std::source_location::current())

#define UV_REQUIRE_STRICTLY_INCREASING(x)                                                \
    ::uv::errors::validateStrictlyIncreasing((x), #x, std::source_location::current())

#define UV_REQUIRE_NON_EMPTY(x)                                                          \
    ::uv::errors::validateNonEmpty((x), #x, std::source_location::current())

#define UV_REQUIRE_SET(x)                                                                \
    ::uv::errors::validateSet((x), #x, std::source_location::current())

#define UV_REQUIRE_SAME_SIZE(a, b)                                                       \
    ::uv::errors::validateSameSize(                                                      \
        (a),                                                                             \
        (b),                                                                             \
        #a " vs " #b,                                                                    \
        std::source_location::current()                                                  \
    )

#define UV_REQUIRE_MIN_SIZE(x, min)                                                      \
    ::uv::errors::validateMinSize((x), (min), #x, std::source_location::current())

#define UV_REQUIRE_DIR_CREATED(ok, path)                                                 \
    ::uv::errors::validateDirCreated((ok), (path), std::source_location::current())

#define UV_REQUIRE_FILE_OPENED(ok, path)                                                 \
    ::uv::errors::validateFileOpened((ok), (path), std::source_location::current())
