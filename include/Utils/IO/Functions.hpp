// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

#include <concepts>
#include <string>

namespace uv::utils
{
/**
 * @brief Pretty-print a 2D matrix with a header row and row labels.
 *
 * Outputs a table-like view via UV_INFO. The first row contains @p title and
 * header entries; each subsequent row prints the corresponding row label
 * followed by the matrix values.
 *
 * @tparam HeaderVec  Iterable container type for the header row.
 * @tparam RowLabels Indexable container type for row labels.
 * @tparam Matrix    2D container type with iterable rows.
 *
 * @param title        Title printed at the top-left of the table.
 * @param header       Header entries printed across the top row.
 * @param rowLabels    Labels printed at the start of each matrix row.
 * @param M            Matrix to print (rows × columns).
 * @param headerPrec   Decimal precision for header values.
 * @param rowLabelPrec Decimal precision for row-label values.
 * @param valuePrec    Decimal precision for matrix values.
 *
 * @note This function assumes compatible dimensions and does not perform
 *       bounds or consistency checks.
 */
template <typename HeaderVec, typename RowLabels, typename Matrix>
void printMatrix(
    std::string_view title,
    const HeaderVec& header,
    const RowLabels& rowLabels,
    const Matrix& M,
    unsigned int headerPrec = 2,
    unsigned int rowLabelPrec = 2,
    unsigned int valuePrec = 5
) noexcept;

/**
 * @brief Pretty-print a 1D vector.
 *
 * Prints all elements of @p v on a single line using fixed-point formatting
 * and the specified decimal precision, and outputs via UV_INFO.
 *
 * @tparam Vector   Iterable container type.
 *
 * @param v         Vector to print.
 * @param valuePrec Decimal precision for values.
 */
template <typename Vector>
void printVector(const Vector& v, unsigned int valuePrec = 5) noexcept;

} // namespace uv::utils

#include "Functions.inl"
