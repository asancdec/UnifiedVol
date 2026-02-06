// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include <Base/Errors/Errors.hpp>

#include <cmath>
#include <cstddef>
#include <format>

namespace uv::errors
{
template <std::floating_point T>
void validateFinite(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        if (!std::isfinite(xs[i])) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] is {}", what, i, xs[i]),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateFinite(T x, std::string_view what, std::source_location loc)
{
    if (!std::isfinite(x)) [[unlikely]]
    {
        raise(ErrorCode::InvalidArgument, std::format("{} is {}", what, x), loc);
    }
}

template <std::floating_point T>
void validateNonNegative(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v = xs[i];

        if (v < 0.0) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be >= 0 (value = {})", what, i, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateNonNegative(T x, std::string_view what, std::source_location loc)
{
    if (x < 0.0) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be >= 0 (value = {})", what, x),
            loc
        );
    }
}

template <std::floating_point T>
void validatePositive(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v = xs[i];

        if (v <= 0.0) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be > 0 (value = {})", what, i, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validatePositive(T x, std::string_view what, std::source_location loc)
{
    if (x <= 0.0) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be > 0 (value = {})", what, x),
            loc
        );
    }
}

template <std::floating_point T>
void validateEqualOrGreater(
    std::span<const T> threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v{xs[i]};
        const T t{threshold[i]};

        if (v < t) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be >= {} (value = {})", what, i, t, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateEqualOrGreater(
    T threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v = xs[i];

        if (v < threshold) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be >= {} (value = {})", what, i, threshold, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateEqualOrGreater(
    T threshold,
    T x,
    std::string_view what,
    std::source_location loc
)
{
    if (x < threshold) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be >= {} (value = {})", what, threshold, x),
            loc
        );
    }
}

template <std::floating_point T>
void validateGreater(
    T threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v = xs[i];

        if (v <= threshold) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be > {} (value = {})", what, i, threshold, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateGreater(T threshold, T x, std::string_view what, std::source_location loc)
{
    if (x <= threshold) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be > {} (value = {})", what, threshold, x),
            loc
        );
    }
}

template <std::floating_point T>
void validateStrictlyIncreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    const std::size_t n{xs.size()};

    T prev{xs.front()};

    for (std::size_t i = 1; i < n; ++i)
    {
        const T curr{xs[i]};

        if (curr <= prev) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format(
                    "{} not strictly increasing at {}: {} <= {}",
                    what,
                    i,
                    curr,
                    prev
                ),
                loc
            );
        }

        prev = curr;
    }
}

template <typename T>
void validateNonEmpty(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    if (xs.empty()) [[unlikely]]
    {
        raise(ErrorCode::InvalidArgument, std::format("{} must be non-empty", what), loc);
    }
}

template <typename T>
void validateSet(
    const std::optional<T>& x,
    std::string_view what,
    std::source_location loc
)
{
    if (!x.has_value()) [[unlikely]]
    {
        raise(ErrorCode::InvalidArgument, std::format("{} must be set", what), loc);
    }
}

} // namespace uv::errors
