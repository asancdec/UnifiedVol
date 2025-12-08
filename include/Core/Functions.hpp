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

#include "Utils/Types.hpp"

#include <cstddef>

namespace uv::core
{
    /**
     * @brief Transpose a dense matrix stored as a vector of row vectors.
     *
     * @tparam Real        Element type.
     * @param  input    Matrix with shape [numRows][numCols].
     * @return          Transposed matrix with shape [numCols][numRows].
     */
	Matrix<Real> transposeMatrix(const Matrix<Real>& input);

    /**
     * @brief Generate an evenly spaced 1D grid of Real values.
     *
     * @param bound1  Lower bound of the grid (inclusive).
     * @param bound2  Upper bound of the grid (inclusive).
     * @param steps   Number of grid points (must be >= 1).
     *
     * @return A vector of Real values: [bound1, ..., bound2].
     *
     * @note If `steps <= 1`, the function returns { bound1 }.
    */
    Vector<Real> generateGrid(const Real bound1,
        const Real bound2,
        const size_t steps
    ) noexcept;

} // namespace uv::core