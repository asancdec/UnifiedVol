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

#pragma once

#include <concepts>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <source_location>
#include <span>
#include <string_view>
#include <type_traits>

namespace uv::errors
{
namespace detail
{

template <typename R>
concept ContiguousFloatRange =
    std::ranges::contiguous_range<R> && std::ranges::sized_range<R> &&
    std::floating_point<std::remove_cvref_t<std::ranges::range_value_t<R>>>;

template <typename R>
using RangeValue = std::remove_cvref_t<std::ranges::range_value_t<R>>;
} // namespace detail

template <std::floating_point T>
void validateFinite(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R>
void validateFinite(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateFinite(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateNonNegative(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R>
void validateNonNegative(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateNonNegative(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validatePositive(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validatePositive(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrLess(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrLess(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrLess(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateLess(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateLess(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    std::span<const T> xs,
    std::span<const T> threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateGreater(
    std::span<const T> xs,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateGreater(
    T x,
    T threshold,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateStrictlyIncreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R>
void validateStrictlyIncreasing(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateNonEmpty(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <detail::ContiguousFloatRange R>
void validateNonEmpty(
    const R& xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateNonNull(
    const T* x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename Ptr>
requires requires(const Ptr& p) { p.get(); }
void validateNonNull(
    const Ptr& p,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateSet(
    const std::optional<T>& x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename A>
requires requires(const A& a) { a.size(); }
void validateSameSize(
    const A& a,
    std::size_t b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename B>
requires requires(const B& b) { b.size(); }
void validateSameSize(
    std::size_t a,
    const B& b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void validateSameSize(
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
void validateSameSize(
    const A& a,
    const B& b,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateMinSize(
    std::span<const T> x,
    std::size_t minSize,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void validateState(
    bool ok,
    std::string_view message,
    std::source_location loc = std::source_location::current()
);

void validateFileOpened(
    bool ok,
    const std::filesystem::path& file,
    std::source_location loc = std::source_location::current()
);

void validateDirCreated(
    bool ok,
    const std::filesystem::path& dir,
    std::source_location loc = std::source_location::current()
);

} // namespace uv::errors

#include "Base/Errors/Detail/Validate.inl"
