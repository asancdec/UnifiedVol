// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <type_traits>

namespace uv::errors::validate::detail
{

template <typename R>
concept ContiguousFloatRange =
    std::ranges::contiguous_range<R> && std::ranges::sized_range<R> &&
    std::floating_point<std::remove_cvref_t<std::ranges::range_value_t<R>>>;

template <typename R> using RangeValue =
    std::remove_cvref_t<std::ranges::range_value_t<R>>;

template <typename T>
concept PointerLike = requires(const T& value) { value.get(); };

template <typename T>
concept Sized = requires(const T& value) {
    { std::size(value) } -> std::convertible_to<std::size_t>;
};

} // namespace uv::errors::validate::detail
