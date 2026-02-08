// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 �lvaro S�nchez de Carlos
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

#include "Base/Errors/Errors.hpp"

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

template <detail::ContiguousFloatRange R>
void validateFinite(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateFinite(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
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

template <detail::ContiguousFloatRange R>
void validateNonNegative(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateNonNegative(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
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
void validateEqualOrLess(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v{xs[i]};
        const T t{threshold[i]};

        if (v > t) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be <= {} (value = {})", what, i, t, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateEqualOrLess(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v{xs[i]};

        if (v > threshold) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be <= {} (value = {})", what, i, threshold, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateEqualOrLess(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc
)
{
    if (x > threshold) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be <= {} (value = {})", what, threshold, x),
            loc
        );
    }
}

template <std::floating_point T>
void validateLess(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc
)
{
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const T v{xs[i]};

        if (v >= threshold) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format("{}[{}] must be < {} (value = {})", what, i, threshold, v),
                loc
            );
        }
    }
}

template <std::floating_point T>
void validateLess(T x, T threshold, std::string_view what, std::source_location loc)
{
    if (x >= threshold) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be < {} (value = {})", what, threshold, x),
            loc
        );
    }
}

template <std::floating_point T>
void validateEqualOrGreater(
    std::span<const T> xs,
    std::span<const T> threshold,
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
    std::span<const T> xs,
    T threshold,
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
    T x,
    T threshold,
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
    std::span<const T> xs,
    T threshold,
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
void validateGreater(T x, T threshold, std::string_view what, std::source_location loc)
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

template <detail::ContiguousFloatRange R>
void validateStrictlyIncreasing(
    const R& xs,
    std::string_view what,
    std::source_location loc
)
{
    using T = detail::RangeValue<R>;
    validateStrictlyIncreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
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
        raise(ErrorCode::InvalidState, std::format("{} must be non-empty", what), loc);
    }
}

template <detail::ContiguousFloatRange R>
void validateNonEmpty(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateNonEmpty(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <typename T>
void validateNonNull(const T* x, std::string_view what, std::source_location loc)
{
    if (x == nullptr) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, std::format("{} must not be null", what), loc);
    }
}

template <typename Ptr>
requires requires(const Ptr& p) { p.get(); }
void validateNonNull(const Ptr& p, std::string_view what, std::source_location loc)
{

    validateNonNull(p.get(), what, loc);
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
        raise(ErrorCode::InvalidState, std::format("{} must be set", what), loc);
    }
}

template <typename A>
requires requires(const A& a) { a.size(); }
void validateSameSize(
    const A& a,
    std::size_t b,
    std::string_view what,
    std::source_location loc
)
{
    if (static_cast<std::size_t>(a.size()) != b) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format(
                "{} size mismatch: {} != {}",
                what,
                static_cast<std::size_t>(a.size()),
                b
            ),
            loc
        );
    }
}

template <typename B>
requires requires(const B& b) { b.size(); }
void validateSameSize(
    std::size_t a,
    const B& b,
    std::string_view what,
    std::source_location loc
)
{
    if (a != static_cast<std::size_t>(b.size())) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format(
                "{} size mismatch: {} != {}",
                what,
                a,
                static_cast<std::size_t>(b.size())
            ),
            loc
        );
    }
}

template <typename A, typename B>
requires requires(const A& a, const B& b) {
    a.size();
    b.size();
}
void validateSameSize(
    const A& a,
    const B& b,
    std::string_view what,
    std::source_location loc
)
{
    if (a.size() != b.size()) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} size mismatch: {} != {}", what, a.size(), b.size()),
            loc
        );
    }
}
template <typename T>
void validateMinSize(
    std::span<const T> x,
    std::size_t minSize,
    std::string_view what,
    std::source_location loc
)
{
    if (x.size() < minSize) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format(
                "{} has size {}, but minimum required size is {}",
                what,
                x.size(),
                minSize
            ),
            loc
        );
    }
}

} // namespace uv::errors
