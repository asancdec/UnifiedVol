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
     * @tparam Real     Element type.
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
    Vector<Real> generateGrid(Real bound1,
        Real bound2,
        size_t steps
    ) noexcept;

    /**
     * @brief Compute first forward differences of a vector.
     *
     * This function returns a vector containing the element-wise
     * forward differences of the input vector `v`, defined as
     * @f$ d[i] = v[i+1] - v[i] @f$ for all valid indices.
     *
     * ### Example
     * Input:  `[v0, v1, v2, v3]`
     * Output: `[v1 - v0, v2 - v1, v3 - v2]`
     *
     * @param v Input vector.
     * @return A vector of size `v.size() - 1` containing forward differences.
     *     
     */
    Vector<Real> diff(const Vector<Real>& v);

    /**
     * @brief Return a new vector with every element multiplied by a scalar.
     *
     * @param v Input vector (unchanged).
     * @param x Scalar multiplier.
     * @return A vector where each element equals `v[i] * x`.
     */
    Vector<Real> multiply(const Vector<Real>& v,
        Real x) noexcept;

    /**
     * @brief Compute the element-wise reciprocal of a vector.
     *
     * @param v Input vector.
     * @return A new vector where each entry equals <tt>1 / v[i]</tt>.
     *
     * @note No range checking is performed; division by zero will result
     *       in undefined behavior.
     */
    Vector<Real> reciprocal(const Vector<Real>& v) noexcept;

    /**
     * @brief Add a scalar value to every element of a matrix.
     *
     * @param A Input matrix.
     * @param x Scalar value to add to each matrix element.
     * @return A matrix where each entry equals <tt>A[i][j] + x</tt>.
     */
    Matrix<Real> add(const Matrix<Real>& A,
        Real x) noexcept;

    /**
     * @brief Compute the element-wise (Hadamard) product of two vectors.
     *
     * @param a First input vector.
     * @param b Second input vector.
     * @return A vector where each element equals <tt>a[i] * b[i]</tt>.
     */
    Vector<Real> hadamard(const Vector<Real>& a, const Vector<Real>& b);

    /**
     * @brief Compute the element-wise (Hadamard) product of a matrix and a vector.
     *
     * @param A Input matrix of size `A.size() x A[0].size()`.
     * @param b Input vector whose size must match the number of columns of @p A.
     * @return A new matrix where each row is the Hadamard product of a row of @p A and @p b.
     */
    Matrix<Real> hadamard(const Matrix<Real>& A,
        const Vector<Real>& b);

    /**
     * @brief Compute the element-wise (Hadamard) product of two matrices.
     *
     * @return A matrix of the same size as `A` and `B`, where each element is the
     *         Hadamard (element-wise) product of the corresponding elements of `A` and `B`.
     */
    Matrix<Real> hadamard(const Matrix<Real>& A,
        const Matrix<Real>& B);

    /**
     * @brief Return the minimum element of a vector.
     *
     * Computes and returns the smallest value in `x` using `std::min_element`.
     *
     * @tparam T Value type. Must be comparable with `operator<`.
     *
     * @param x Input vector (must be non-empty).
     *
     * @return Minimum value in `x`.
     *
     * @throws ErrorCode::InvalidArgument if `x` is empty.
     */
    template<typename T>
    T minValue(const Vector<T>& x);

    /**
     * @brief Return the maximum element of a vector.
     *
     * Computes and returns the largest value in `x` using `std::max_element`.
     *
     * @tparam T Value type. Must be comparable with `operator<`.
     *
     * @param x Input vector (must be non-empty).
     *
     * @return Maximum value in `x`.
     *
     * @throws ErrorCode::InvalidArgument if `x` is empty.
     */
    template<typename T>
    T maxValue(const Vector<T>& x);

    /**
     * @brief Convert a vector to another value type.
     *
     * Creates and returns a new vector by converting each element of the input
     * vector `x` to type `To` using `static_cast`.
     *
     * The input vector may be empty; in this case, an empty vector is returned.
     * No shape or validity checks are performed beyond element-wise conversion.
     *
     * @tparam To   Target value type.
     * @tparam From Source value type.
     *
     * @param x Input vector.
     *
     * @return Vector whose elements are `static_cast<To>(x[i])`.
     */
    template <typename To, typename From>
    Vector<To> convertVector(const Vector<From>& x) noexcept;

    /**
     * @brief Convert a matrix to another value type.
     *
     * Creates and returns a new matrix by converting each element of the input
     * matrix `A` to type `To` using `static_cast`. The row structure of the matrix
     * is preserved.
     *
     * The input matrix may be empty, and rows may be empty; these cases are handled
     * naturally and result in an empty matrix or empty rows in the output.
     *
     * @tparam To   Target value type.
     * @tparam From Source value type.
     *
     * @param A Input matrix.
     *
     * @return Matrix whose elements are `static_cast<To>(A[i][j])`.
     */
    template <typename To, typename From>
    Matrix<To> convertMatrix(const Matrix<From>& A) noexcept;


} // namespace uv::core

#include "Functions.inl"