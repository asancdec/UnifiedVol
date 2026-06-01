// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"

#include <cmath>
#include <cstddef>
#include <format>

namespace uv::errors
{
namespace detail
{
template <std::floating_point T> void validateStrictOrder(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc,
    bool increasing
)
{
    const std::size_t n{xs.size()};

    if (n < 2)
        return;

    T prev{xs.front()};

    for (std::size_t i = 1; i < n; ++i)
    {
        const T curr{xs[i]};
        const bool ok{increasing ? curr > prev : curr < prev};

        if (!ok) [[unlikely]]
        {
            raise(
                ErrorCode::InvalidArgument,
                std::format(
                    "{} not strictly {} at {}: {} {} {}",
                    what,
                    increasing ? "increasing" : "decreasing",
                    i,
                    curr,
                    increasing ? "<=" : ">=",
                    prev
                ),
                loc
            );
        }

        prev = curr;
    }
}
} // namespace detail

template <std::floating_point T> void
validateFinite(std::span<const T> xs, std::string_view what, std::source_location loc)
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

template <std::floating_point T> void validateNonNegative(
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

template <std::floating_point T> void
validatePositive(std::span<const T> xs, std::string_view what, std::source_location loc)
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

template <class A, class B>
requires std::equality_comparable_with<A, B>
void validateEqual(
    const A& a,
    const B& b,
    std::string_view what,
    std::source_location loc
)
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
void validateClose(T a, T b, T tol, std::string_view what, std::source_location loc)
{
    validateNonNegative(tol, "tol", loc);

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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <class T>
requires std::totally_ordered<T>
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

template <std::floating_point T> void validateStrictlyIncreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    detail::validateStrictOrder(xs, what, loc, true);
}

template <detail::ContiguousFloatRange R> void
validateStrictlyIncreasing(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateStrictlyIncreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T> void validateStrictlyDecreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    detail::validateStrictOrder(xs, what, loc, false);
}

template <detail::ContiguousFloatRange R> void
validateStrictlyDecreasing(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateStrictlyDecreasing(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <std::floating_point T> void validateStrictlyMonotonic(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc
)
{
    const std::size_t n{xs.size()};

    if (n < 2)
        return;

    if (xs[1] > xs[0])
    {
        validateStrictlyIncreasing(xs, what, loc);
        return;
    }

    if (xs[1] < xs[0])
    {
        validateStrictlyDecreasing(xs, what, loc);
        return;
    }

    raise(
        ErrorCode::InvalidArgument,
        std::format("{} not strictly monotonic at 1: {} == {}", what, xs[1], xs[0]),
        loc
    );
}

template <detail::ContiguousFloatRange R> void
validateStrictlyMonotonic(const R& xs, std::string_view what, std::source_location loc)
{
    using T = detail::RangeValue<R>;
    validateStrictlyMonotonic(
        std::span<const T>(std::ranges::data(xs), std::ranges::size(xs)),
        what,
        loc
    );
}

template <typename T> void
validateNonEmpty(std::span<const T> xs, std::string_view what, std::source_location loc)
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

template <typename T> void
validateSet(const std::optional<T>& x, std::string_view what, std::source_location loc)
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
template <typename T> void validateMinSize(
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
