// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/Generate.hpp"
#include "Core/MarketData.hpp"
#include "Core/MarketState.hpp"
#include "IO/CSV/Read.hpp"

#include <concepts>
#include <filesystem>
#include <string>
#include <tuple>

namespace uv::io::load
{
template <std::floating_point T> core::MarketState<T> marketState(
    const std::filesystem::path& path,
    const core::MarketData<T>& marketData,
    const csv::Options& opt = {}
);

namespace detail
{
template <std::floating_point T> std::tuple<Vector<T>, Vector<T>, core::Matrix<T>>
readLabeledMatrixCsv(const std::string& filename, const csv::Options& opt = {});
}
} // namespace uv::io::load

#include "IO/Detail/Load.inl"