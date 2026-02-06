// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
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

#include <Base/Alias.hpp>
#include <Core/Generate.hpp>
#include <Core/MarketData.hpp>
#include <Core/MarketState.hpp>
#include <IO/CSV/Read.hpp>

#include <concepts>
#include <filesystem>
#include <string>
#include <tuple>

namespace uv::io::load
{
template <std::floating_point T>
core::MarketState<T> marketState(
    const std::filesystem::path& path,
    const core::MarketData<T>& marketData,
    const csv::Options& opt = {}
);

namespace detail
{
template <std::floating_point T>
std::tuple<Vector<T>, Vector<T>, core::Matrix<T>>
readLabeledMatrixCsv(const std::string& filename, const csv::Options& opt = {});
}
} // namespace uv::io::load

#include <IO/Detail/Load.inl>