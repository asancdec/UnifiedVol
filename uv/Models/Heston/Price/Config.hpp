// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <cstddef>

namespace uv::models::heston::price
{
template <std::floating_point T> struct Config
{
    T alphaItm{T{-2}};
    T alphaOtm{T{2}};
};

inline constexpr std::size_t defaultNodes{500};

} // namespace uv::models::heston::price
