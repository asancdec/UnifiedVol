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

#include "Base/Macros/Require.hpp"

#include <fstream>
#include <vector>

namespace uv::io::load
{

template <std::floating_point T>
core::MarketState<T> marketState(
    const std::filesystem::path& path,
    const core::MarketData<T>& marketData,
    const csv::Options& opt
)
{
    auto [maturities, moneyness, vol] =
        detail::readLabeledMatrixCsv<T>(path.string(), opt);

    return core::generateMarketState<T>(marketData, maturities, moneyness, vol);
}
} // namespace uv::io::load

namespace uv::io::load::detail
{
template <std::floating_point T>
std::tuple<Vector<T>, Vector<T>, core::Matrix<T>>
readLabeledMatrixCsv(const std::string& filename, const csv::Options& opt)
{
    std::ifstream file(filename);

    UV_REQUIRE_FILE_OPENED(file.is_open(), filename);

    auto dense = csv::readLabeledDenseOrThrow<T, Vector>(file, filename, opt);

    core::Matrix<T> mat(dense.rows, dense.cols);

    for (std::size_t i = 0; i < dense.rows; ++i)
        for (std::size_t j = 0; j < dense.cols; ++j)
            mat[i][j] = dense.values[i * dense.cols + j];

    return {std::move(dense.rowLabels), std::move(dense.colLabels), std::move(mat)};
}
} // namespace uv::io::load::detail
