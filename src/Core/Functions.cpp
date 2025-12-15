// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.cpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include "Core/Functions.hpp"
#include "Utils/Aux/Errors.hpp"  
#include "Utils/Types.hpp"

#include <cstddef>
#include <string> 
#include <vector>

namespace uv::core
{
    Matrix<Real> transposeMatrix(const Matrix<Real>& input)
	{
        // ---------- Check matching dimensions ----------

        const std::size_t numRows{ input.size() };
        const std::size_t numCols{ input[0].size() };

        // Throw if all rows are not the same size
        for (std::size_t i = 1; i < numRows; ++i)
        {
            const std::size_t actual{ input[i].size() };
            const std::size_t expected{ numCols };

            UV_REQUIRE(
                actual == expected,
                ErrorCode::InvalidArgument,
                "transposeVector: inconsistent row length - row " +
                std::to_string(i) + " has " + std::to_string(actual) +
                " elements, expected " + std::to_string(expected)
            );
        }

        // ---------- Transpose matrix ----------
        
        std::vector<Vector<Real>> result(numCols, Vector<Real>(numRows));

        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                result[j][i] = input[i][j];
            }
        }

        return result;
	}

    Vector<Real> generateGrid(const Real bound1,
        const Real bound2,
        const size_t steps
    ) noexcept
    {
        // ---------- Allocate memory ----------

        Vector<Real> grid;
        grid.reserve(steps);

        // ---------- Handle special case ----------

        if (steps <= 1)
        {
            grid.push_back(static_cast<Real>(bound1));
            return grid;
        }

        // ---------- Generate grid ----------

        const Real dx{ (bound2 - bound1) / static_cast<Real>(steps - 1) };

        for (std::size_t i = 0; i < steps; ++i)
        {
            grid.push_back(static_cast<Real>(bound1) + dx * static_cast<Real>(i));
        }

        return grid;
    }

    Vector<Real> diff(const Vector<Real>& v)
    {
        // ---------- Check dimensions ----------

        const std::size_t n{ v.size() };

        UV_REQUIRE(
            n >= 2,
            ErrorCode::InvalidArgument,
            "diff: input vector must contain at least 2 elements (got " +
            std::to_string(v.size()) + ")"
        );

        // ---------- Allocate ----------

        // One element less
        Vector<Real> d(n - 1);

        // ---------- Calculate differences ----------

        for (std::size_t i = 0; i + 1 < n; ++i)
        {
            d[i] = v[i + 1] - v[i];
        }
        return d;
    }

    Vector<Real> multiply(const Vector<Real>& v, 
        const Real x) noexcept
    {
        std::size_t vSize{v.size()};

        Vector<Real> result(vSize);

        for (std::size_t i = 0; i < vSize; ++i)
        {
            result[i] = v[i] * x;
        }

        return result;
    }

    Vector<Real> reciprocal(const Vector<Real>& v) noexcept
    {
        std::size_t vSize{ v.size() };

        Vector<Real> result(vSize);

        for (std::size_t i = 0; i < vSize; ++i)
        {
            result[i] = Real(1.0) / v[i];
        }

        return result;
    }


    Matrix<Real> add(const Matrix<Real>& A,
        const Real x) noexcept
    {
        const std::size_t rows{ A.size() };

        Matrix<Real> result(rows);

        for (std::size_t i = 0; i < rows; ++i)
        {
            const std::size_t cols{ A[i].size() };
            result[i].resize(cols);

            for (std::size_t j = 0; j < cols; ++j)
            {
                result[i][j] = A[i][j] + x;
            }
        }

        return result;
    }

    Vector<Real> hadamard(const Vector<Real>& a, 
        const Vector<Real>& b)
    {
        // ---------- Check dimensions ----------

        std::size_t aSize{ a.size() };
        std::size_t bSize{ b.size() };

        UV_REQUIRE(
            aSize == bSize,
            ErrorCode::InvalidArgument,
            "hadamard: vectors must have same size (got " +
            std::to_string(aSize) + " and " + std::to_string(bSize) + ")"
        );

        // ---------- Element-wise multiplication ----------

        Vector<Real> c(aSize);

        for (std::size_t i = 0; i < aSize; ++i)
        {
            c[i] = a[i] * b[i];
        }

        return c;
    }


    Matrix<Real> hadamard(const Matrix<Real>& A,
        const Vector<Real>& b) 
    {
        // ---------- Check dimensions ----------

        const std::size_t rows{ A.size() };
        const std::size_t cols{ A[0].size() };

        UV_REQUIRE(
            b.size() == cols,
            ErrorCode::InvalidArgument,
            "hadamard: vector size must equal number of columns"
        );

        // ---------- Row-wise multiplication ----------

        Matrix<Real> C(rows, Vector<Real>(cols));

        for (std::size_t i = 0; i < rows; ++i)
        {

            C[i] = hadamard(A[i], b);
        }

        return C;
    }

    Matrix<Real> hadamard(const Matrix<Real>& A,
        const Matrix<Real>& B)
    {
        // ---------- Check dimensions ----------

        const std::size_t rowsA{ A.size() };
        const std::size_t rowsB{ B.size() };

        UV_REQUIRE(
            rowsA == rowsB,
            ErrorCode::InvalidArgument,
            "hadamard: matrices must have same number of rows (got " +
            std::to_string(rowsA) + " and " + std::to_string(rowsB) + ")"
        );

        const std::size_t colsA{ A[0].size() };
        const std::size_t colsB{ B[0].size() };

        UV_REQUIRE(
            colsA == colsB,
            ErrorCode::InvalidArgument,
            "hadamard: matrices must have same number of columns (got " +
            std::to_string(colsA) + " and " + std::to_string(colsB) + ")"
        );

        // ---------- Element-wise multiplication ----------

        Matrix<Real> C(rowsA, Vector<Real>(colsA));

        for (std::size_t i = 0; i < rowsA; ++i)
        {
            for (std::size_t j = 0; j < colsA; ++j)
            {
                C[i][j] = A[i][j] * B[i][j];
            }
        }

        return C;
    }

} // namespace uv::core