// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
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

} // namespace uv::errors::validate::detail
