// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/MarketData.hpp"
#include "Core/MarketState.hpp"
#include "IO/CSV/Detail/Read.hpp"

#include <concepts>
#include <filesystem>

namespace uv::io::csv::load
{
using Options = detail::Options;

template <std::floating_point T> core::MarketState<T> marketState(
    const std::filesystem::path& path,
    const core::MarketData<T>& marketData,
    const Options& opt = {}
);
} // namespace uv::io::csv::load

#include "IO/CSV/Detail/Load.inl"
