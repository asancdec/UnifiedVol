// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"

#include <cmath>
#include <cstddef>
#include <format>

namespace uv::errors::validate
{
namespace detail
{
template <std::floating_point T> void strictOrder(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc,
    bool increasing,
    bool strict
)
{
    const std::size_t n{xs.size()};

    if (n < 2)
        return;

    T prev{xs.front()};

    for (std::size_t i = 1; i < n; ++i)
    {
        const T curr{xs[i]};
        const bool ok{
            strict ? (increasing ? curr > prev : curr < prev)
                   : (increasing ? curr >= prev : curr <= prev)
        };

        if (!ok) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format(
                    "{} not {}{} at {}: {} {} {}",
                    what,
                    strict ? "strictly " : "",
                    increasing ? "increasing" : "decreasing",
                    i,
                    curr,
                    increasing ? (strict ? "<=" : "<") : (strict ? ">=" : ">"),
                    prev
                ),
                loc
            );
        }

        prev = curr;
    }
}
} // namespace detail

template <std::floating_point T>
void finite(std::span<const T> xs, std::string_view what, std::source_location loc)
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
void finite(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    finite(std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)), what, loc);
}

template <std::floating_point T>
void finite(T x, std::string_view what, std::source_location loc)
{
    if (!std::isfinite(x)) [[unlikely]]
    {
        raise(ErrorCode::InvalidArgument, std::format("{} is {}", what, x), loc);
    }
}

template <std::floating_point T>
void nonNegative(std::span<const T> xs, std::string_view what, std::source_location loc)
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
void nonNegative(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    nonNegative(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T>
void nonNegative(T x, std::string_view what, std::source_location loc)
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
void positive(std::span<const T> xs, std::string_view what, std::source_location loc)
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

template <detail::ContiguousFloatRange R>
void positive(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    positive(std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)), what, loc);
}

template <std::floating_point T>
void positive(T x, std::string_view what, std::source_location loc)
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

template <class A, class B>
requires std::equality_comparable_with<A, B>
void equal(const A& a, const B& b, std::string_view what, std::source_location loc)
{
    if (a != b) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be equal ({} != {})", what, a, b),
            loc
        );
    }
}

template <std::floating_point T>
void close(T a, T b, T tol, std::string_view what, std::source_location loc)
{
    nonNegative(tol, "tol", loc);

    const T diff{std::abs(a - b)};

    if (diff > tol) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} must be close: |{} - {}| = {} > {}", what, a, b, diff, tol),
            loc
        );
    }
}

template <class T>
requires std::totally_ordered<T>
void equalOrLess(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc
)
{
    sameSize(xs, threshold, what, loc);

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

template <class T>
requires std::totally_ordered<T>
void equalOrLess(
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

template <class T>
requires std::totally_ordered<T>
void equalOrLess(T x, T threshold, std::string_view what, std::source_location loc)
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

template <class T>
requires std::totally_ordered<T>
void less(
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

template <class T>
requires std::totally_ordered<T>
void less(T x, T threshold, std::string_view what, std::source_location loc)
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

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc
)
{
    sameSize(xs, threshold, what, loc);

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

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(
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

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(T x, T threshold, std::string_view what, std::source_location loc)
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

template <class T>
requires std::totally_ordered<T>
void greater(
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

template <class T>
requires std::totally_ordered<T>
void greater(T x, T threshold, std::string_view what, std::source_location loc)
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

template <std::floating_point T> void
strictlyIncreasing(std::span<const T> xs, std::string_view what, std::source_location loc)
{
    detail::strictOrder(xs, what, loc, true, true);
}

template <detail::ContiguousFloatRange R>
void strictlyIncreasing(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    strictlyIncreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T>
void nonDecreasing(std::span<const T> xs, std::string_view what, std::source_location loc)
{
    detail::strictOrder(xs, what, loc, true, false);
}

template <detail::ContiguousFloatRange R>
void nonDecreasing(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    nonDecreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T> void
strictlyDecreasing(std::span<const T> xs, std::string_view what, std::source_location loc)
{
    detail::strictOrder(xs, what, loc, false, true);
}

template <detail::ContiguousFloatRange R>
void strictlyDecreasing(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    strictlyDecreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T> void
strictlyMonotonic(std::span<const T> xs, std::string_view what, std::source_location loc)
{
    const std::size_t n{xs.size()};

    if (n < 2)
        return;

    if (xs[1] > xs[0])
    {
        strictlyIncreasing(xs, what, loc);
        return;
    }

    if (xs[1] < xs[0])
    {
        strictlyDecreasing(xs, what, loc);
        return;
    }

    raise(
        ErrorCode::InvalidArgument,
        std::format("{} not strictly monotonic at 1: {} == {}", what, xs[1], xs[0]),
        loc
    );
}

template <detail::ContiguousFloatRange R>
void strictlyMonotonic(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    strictlyMonotonic(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <typename T>
void nonEmpty(std::span<const T> xs, std::string_view what, std::source_location loc)
{
    if (xs.empty()) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, std::format("{} must be non-empty", what), loc);
    }
}

template <detail::ContiguousFloatRange R>
void nonEmpty(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    nonEmpty(std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)), what, loc);
}

template <typename T>
void nonNull(const T* x, std::string_view what, std::source_location loc)
{
    if (x == nullptr) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, std::format("{} must not be null", what), loc);
    }
}

template <typename Ptr>
requires requires(const Ptr& p) { p.get(); }
void nonNull(const Ptr& p, std::string_view what, std::source_location loc)
{

    nonNull(p.get(), what, loc);
}

template <typename T>
void set(const std::optional<T>& x, std::string_view what, std::source_location loc)
{
    if (!x.has_value()) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, std::format("{} must be set", what), loc);
    }
}

template <typename A>
requires requires(const A& a) { a.size(); }
void sameSize(const A& a, std::size_t b, std::string_view what, std::source_location loc)
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
void sameSize(std::size_t a, const B& b, std::string_view what, std::source_location loc)
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
void sameSize(const A& a, const B& b, std::string_view what, std::source_location loc)
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
template <typename T> void minSize(
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

} // namespace uv::errors::validate
