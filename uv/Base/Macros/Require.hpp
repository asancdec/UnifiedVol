// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Errors/Errors.hpp"   // IWYU pragma: keep
#include "Base/Errors/Validate.hpp" // IWYU pragma: keep

#include <source_location> // IWYU pragma: keep

#define REQUIRE(cond, code, message)                                                     \
    do                                                                                   \
    {                                                                                    \
        if (!(cond))                                                                     \
        {                                                                                \
            ::uv::errors::raise((code), (message), std::source_location::current());     \
        }                                                                                \
    } while (0)

#define REQUIRE_FINITE(x)                                                                \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::finite((x), #x, std::source_location::current());        \
    } while (0)

#define REQUIRE_NON_NEGATIVE(x)                                                          \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::nonNegative((x), #x, std::source_location::current());   \
    } while (0)

#define REQUIRE_POSITIVE(x)                                                              \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::positive((x), #x, std::source_location::current());      \
    } while (0)

#define REQUIRE_EQUAL(a, b)                                                              \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::equal(                                                   \
            (a),                                                                         \
            (b),                                                                         \
            #a " vs " #b,                                                                \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_CLOSE(a, b, ...)                                                         \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::close(                                                   \
            (a),                                                                         \
            (b),                                                                         \
            __VA_ARGS__ __VA_OPT__(, ) #a " vs " #b,                                     \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_EQUAL_OR_LESS(x, threshold)                                              \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::equalOrLess(                                             \
            (x),                                                                         \
            (threshold),                                                                 \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_LESS(x, threshold)                                                       \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::less(                                                    \
            (x),                                                                         \
            (threshold),                                                                 \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_EQUAL_OR_GREATER(x, threshold)                                           \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::equalOrGreater(                                          \
            (x),                                                                         \
            (threshold),                                                                 \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_GREATER(x, threshold)                                                    \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::greater(                                                 \
            (x),                                                                         \
            (threshold),                                                                 \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_STRICTLY_INCREASING(x)                                                   \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::strictlyIncreasing(                                      \
            (x),                                                                         \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_NON_DECREASING(x)                                                        \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::nonDecreasing((x), #x, std::source_location::current()); \
    } while (0)

#define REQUIRE_STRICTLY_DECREASING(x)                                                   \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::strictlyDecreasing(                                      \
            (x),                                                                         \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_STRICTLY_MONOTONIC(x)                                                    \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::strictlyMonotonic(                                       \
            (x),                                                                         \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_NON_EMPTY(x)                                                             \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::nonEmpty((x), #x, std::source_location::current());      \
    } while (0)

#define REQUIRE_NON_NULL(x)                                                              \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::nonNull((x), #x, std::source_location::current());       \
    } while (0)

#define REQUIRE_SET(x)                                                                   \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::set((x), #x, std::source_location::current());           \
    } while (0)

#define REQUIRE_SAME_SIZE(a, b)                                                          \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::sameSize(                                                \
            (a),                                                                         \
            (b),                                                                         \
            #a " vs " #b,                                                                \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_MIN_SIZE(x, minSizeValue)                                                \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::minSize(                                                 \
            (x),                                                                         \
            (minSizeValue),                                                              \
            #x,                                                                          \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_VALID_STATE(ok, message)                                                 \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::state((ok), (message), std::source_location::current()); \
    } while (0)

#define REQUIRE_DIR_CREATED(ok, path)                                                    \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::dirCreated(                                              \
            (ok),                                                                        \
            (path),                                                                      \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)

#define REQUIRE_FILE_OPENED(ok, path)                                                    \
    do                                                                                   \
    {                                                                                    \
        ::uv::errors::validate::fileOpened(                                              \
            (ok),                                                                        \
            (path),                                                                      \
            std::source_location::current()                                              \
        );                                                                               \
    } while (0)
