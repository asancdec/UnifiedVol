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

#include "Core/VolSurface.hpp"

#include <string>

namespace uv::utils
{
	/**
	 * @brief Reads a CSV file containing a volatility surface and returns a VolSurface object.
	 *
	 * The CSV file is expected to have the following structure:
	 * - The first row contains moneyness (K/S).
	 * - The first column of each subsequent row contains tenors.
	 * - The remaining cells contain implied volatilities for each strike-maturity pair.
	 *
	 * @param filename Path to the CSV file.
	 * @return VolSurface Object containing the strikes, tenors, and implied volatilities.
	 */
	core::VolSurface readVolSurface(const std::string& filename, const core::MarketData& mktData);

	/**
	 * @brief Pretty-print a 2D matrix with a header row and row labels.
	 *
	 * Prints a table-like view to the console/logger via UV_INFO. The first line
	 * contains @p title, followed by a header row (e.g., strikes or moneyness).
	 * Each subsequent line prints the corresponding row label (e.g., tenor) and
	 * the values of that matrix row.
	 *
	 * Expected layout:
	 * - Row 0: title + header entries
	 * - Rows i>0: rowLabels[i] + M[i][*]
	 *
	 * @tparam HeaderVec  Iterable container type for the header row.
	 * @tparam RowLabels  Indexable container type for row labels (operator[]).
	 * @tparam Matrix     2D container type where M.size() gives row count and
	 *                    M[i] is iterable over the row values.
	 *
	 * @param title       Title printed at the top-left of the table.
	 * @param header      Header entries printed across the top row.
	 * @param rowLabels   Labels printed at the start of each matrix row.
	 * @param M           Matrix to print (rows x columns).
	 * @param headerPrec  Decimal precision for header values.
	 * @param rowLabelPrec Decimal precision for row-label values.
	 * @param valuePrec   Decimal precision for matrix values.
	 *
	 * @note This function does not validate dimensions. It assumes:
	 *       - rowLabels has at least M.size() elements
	 *       - each M[i] is iterable
	 */
	template <typename HeaderVec, typename RowLabels, typename Matrix>
	void printMatrix(std::string_view title,
		const HeaderVec& header,
		const RowLabels& rowLabels,
		const Matrix& M,
		unsigned int headerPrec = 2,
		unsigned int rowLabelPrec = 2,
		unsigned int valuePrec = 5) noexcept;

} // namespace uv::utils

#include "Functions.inl"

