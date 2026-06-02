// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Errors/Detail/ValidateConcepts.hpp"

#include <concepts>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <source_location>
#include <span>
#include <string_view>

namespace uv::errors::validate
{

template <std::floating_point T> void finite(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void finite(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void finite(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void nonNegative(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void nonNegative(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void nonNegative(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void positive(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void positive(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void positive(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class A, class B>
requires std::equality_comparable_with<A, B>
void equal(
    const A& a,
    const B& b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void close(
    T a,
    T b,
    T tol = T{1e-10},
    std::string_view what = {},
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrLess(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrLess(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrLess(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void less(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void less(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void equalOrGreater(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void greater(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <class T>
requires std::totally_ordered<T>
void greater(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void strictlyIncreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void strictlyIncreasing(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void strictlyDecreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void strictlyDecreasing(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T> void strictlyMonotonic(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void strictlyMonotonic(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T> void nonEmpty(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R> void nonEmpty(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T> void nonNull(
    const T* x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename Ptr>
requires requires(const Ptr& p) { p.get(); }
void nonNull(
    const Ptr& p,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T> void
set(const std::optional<T>& x,
    std::string_view what,
    std::source_location loc = std::source_location::current());

template <typename A>
requires requires(const A& a) { a.size(); }
void sameSize(
    const A& a,
    std::size_t b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename B>
requires requires(const B& b) { b.size(); }
void sameSize(
    std::size_t a,
    const B& b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void sameSize(
    std::size_t a,
    std::size_t b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename A, typename B>
requires requires(const A& a, const B& b) {
    a.size();
    b.size();
}
void sameSize(
    const A& a,
    const B& b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T> void minSize(
    std::span<const T> x,
    std::size_t minSize,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void state(
    bool ok,
    std::string_view message,
    std::source_location loc = std::source_location::current()
);

void fileOpened(
    bool ok,
    const std::filesystem::path& file,
    std::source_location loc = std::source_location::current()
);

void dirCreated(
    bool ok,
    const std::filesystem::path& dir,
    std::source_location loc = std::source_location::current()
);

} // namespace uv::errors::validate

#include "Base/Errors/Detail/Validate.inl"
