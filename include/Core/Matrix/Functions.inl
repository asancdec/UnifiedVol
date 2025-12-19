// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
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

#include "Utils/Aux/Errors.hpp"

#include <cstddef>
#include <span>
#include <cmath>

namespace uv::core
{
    template <std::floating_point T, typename F>
    requires std::invocable<F&, std::size_t, std::size_t>
    Matrix<T> generateIndexed(std::size_t rows,
        std::size_t cols,
        F&& f)
    {
        Matrix<T> out(rows, cols);

        for (std::size_t i = 0; i < rows; ++i)
        {
            std::span<T> row{ out[i] };

            for (std::size_t j = 0; j < cols; ++j)
            {
                row[j] = std::invoke(f, i, j);
            }
        }
        return out;
    }

    template <std::floating_point T, typename F>
    requires std::invocable<F&, std::size_t, std::size_t, T>
    Matrix<T> transformIndexed(const Matrix<T>& m, F&& f)
    {
        return generateIndexed<T>
            (
            m.rows(), 
            m.cols(),
            [&](std::size_t i, std::size_t j)
            {
                return std::invoke(f, i, j, m[i][j]);
            }
        );
    }

    template <std::floating_point T, typename F>
    requires std::invocable<F&, std::size_t, std::size_t, T>
    void transformIndexedInplace(Matrix<T>& m, F&& f)
    {
        for (std::size_t i = 0; i < m.rows(); ++i)
        {
            std::span<T> row{ m[i] };

            for (std::size_t j = 0; j < m.cols(); ++j)
            {
                row[j] = std::invoke(f, i, j, row[j]);
            }  
        }
    }

    template <std::floating_point T>
    Matrix<T> hadamard(const Matrix<T>& lhs,
        const Matrix<T>& rhs)
    {
        Matrix<T> out{ lhs };       
        hadamardInplace(out, rhs);  
        return out;
    }

    template <std::floating_point T>
    void hadamardInplace(Matrix<T>& lhs,
        const Matrix<T>& rhs)
    {
        const std::size_t numRows{ lhs.rows() };
        const std::size_t numCols{ lhs.cols() };

        UV_REQUIRE
        (
            numRows == rhs.rows() && numCols == rhs.cols(),
            ErrorCode::InvalidArgument,
            "hadamardInplace: matrix dimensions must match"
        );

        for (std::size_t i = 0; i < numRows; ++i)
        {
            std::span<T> lhsRow{ lhs[i] };
            std::span<const T> rhsRow {rhs[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                lhsRow[j] *= rhsRow[j];
            }
        }
    }

    template <std::floating_point T>
    Matrix<T> hadamard(
        const Matrix<T>& lhs,
        const Vector<T>& rhs
    )
    {
        Matrix<T> out{ lhs };          
        hadamardInplace(out, rhs);     
        return out;
    }

    template <std::floating_point T>
    void hadamardInplace(Matrix<T>& lhs,
        const Vector<T>& rhs)
    {
        const std::size_t numRows{ lhs.rows() };
        const std::size_t numCols{ lhs.cols() };

        UV_REQUIRE(
            numRows == rhs.size(),
            ErrorCode::InvalidArgument,
            "hadamardInplace: vector size must match number of rows"
        );

        for (std::size_t i = 0; i < numRows; ++i)
        {
            std::span<T> row{ lhs[i] };
            const T scale{ rhs[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                row[j] *= scale;
            }
        }
    }

    template <std::floating_point T>
    Matrix<T> divide(const Matrix<T>& lhs,
        const Matrix<T>& rhs)
    {
        Matrix<T> out{lhs};
        divideInplace(out, rhs);
        return out;
    }

    template <std::floating_point T>
    void divideInplace(Matrix<T>& lhs,
        const Matrix<T>& rhs)
    {
        const std::size_t numRows{ lhs.rows() };
        const std::size_t numCols{ lhs.cols() };

        UV_REQUIRE
        (
            numRows == rhs.rows() && numCols == rhs.cols(),
            ErrorCode::InvalidArgument,
            "divideInplace: matrix dimensions must match"
        );

        for (std::size_t i = 0; i < numRows; ++i)
        {
            std::span<T> lhsRow{ lhs[i] };
            std::span<const T> rhsRow{ rhs[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                lhsRow[j] /= rhsRow[j];
            }
        }
    }

    template<std::floating_point T>
    Matrix<T> reciprocal(const Matrix<T>& m) noexcept
    {
        Matrix<T> out{ m };
        reciprocalInplace(out);
        return out;
    }

    template<std::floating_point T>
    void reciprocalInplace(Matrix<T>& m) noexcept
    {
        const std::size_t numCols{ m.cols() };

        for (std::size_t i = 0; i <  m.rows(); ++i)
        {
            std::span<T> mRow{ m[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                mRow[j] = T{ 1 } / mRow[j];
            }
        }
    }
 
    template <std::floating_point T>
    Matrix<T> square(const Matrix<T>& m) noexcept
    {
        Matrix<T> out{ m };
        squareInplace(out);
        return out;
    }

    template <std::floating_point T>
    void squareInplace(Matrix<T>& m) noexcept
    {
        const std::size_t numCols{ m.cols() };

        for (std::size_t i = 0; i < m.rows(); ++i)
        {
            std::span<T> row{ m[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                row[j] *= row[j];
            }
        }
    }

    template <std::floating_point T>
    Matrix<T> sqrt(const Matrix<T>& m)
    {
        Matrix<T> out{ m };
        sqrtInplace(out);
        return out;
    }

    template <std::floating_point T>
    void sqrtInplace(Matrix<T>& m)
    {
        const std::size_t numRows{ m.rows() };
        const std::size_t numCols{ m.cols() };

        for (std::size_t i = 0; i < numRows; ++i)
        {
            std::span<T> row{ m[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                UV_REQUIRE(
                    row[j] >= T{ 0 },
                    ErrorCode::InvalidArgument,
                    "sqrtInplace: negative matrix entry"
                );

                row[j] = std::sqrt(row[j]);
            }
        }
    }

    template <typename To, typename From>
    Matrix<To> convertMatrix(const Matrix<From>& A) noexcept
    {
        std::size_t numRows{ A.rows() };
        std::size_t numCols{ A.cols() };

        core::Matrix<To> out(numRows, numCols);

        for (std::size_t i = 0; i < numRows; ++i)
        {
            std::span<const From> inRow{ A[i] };
            std::span<To> outRow{ out[i] };

            for (std::size_t j = 0; j < numCols; ++j)
            {
                outRow[j] = static_cast<To>(inRow[j]);
            }
        }
        return out;
    }
} // namespace uv::core