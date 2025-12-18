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

#include "Core/Matrix/Matrix.hpp"
#include "Utils/Types.hpp"

#include <concepts>

namespace uv::core
{
    /**
     * @brief Apply a function to each matrix element using its indices.
     *
     * Returns a new matrix where each element is computed as f(i, j, m[i][j]).
     */
    template <std::floating_point T, typename F>
    Matrix<T> applyIndexed(const Matrix<T>& m, F&& f);

    /**
     * @brief Apply a function to each matrix element in place using its indices.
     *
     * Each element m[i][j] is replaced by f(i, j, m[i][j]).
     */
    template <std::floating_point T, typename F>
    void applyIndexedInplace(Matrix<T>& m, F&& f);

    /**
     * @brief Create a matrix by applying a function over index pairs.
     *
     * Returns a rows by cols matrix with elements f(i, j).
     */
    template <std::floating_point T, typename F>
    Matrix<T> applyIndexed(std::size_t rows,
        std::size_t cols,
        F&& f
    );

    /**
     * @brief Fill a matrix in place by applying a function over index pairs.
     *
     * Each element m[i][j] is replaced by f(i, j).
     */
    template <std::floating_point T, typename F>
    void applyIndexedInplace(Matrix<T>& m,
        std::size_t rows,
        std::size_t cols,
        F&& f
    );

    /**
     * @brief Element wise multiplication of two matrices.
     *
     * Returns a new matrix where each element is lhs[i][j] * rhs[i][j].
     */
    template<std::floating_point T>
    Matrix<T> hadamard(const Matrix<T>& lhs,
        const Matrix<T>& rhs);

    /**
     * @brief In place element wise multiplication with another matrix.
     *
     * Replaces lhs[i][j] with lhs[i][j] * rhs[i][j].
     */
    template<std::floating_point T>
    void hadamardInplace(Matrix<T>& lhs,
        const Matrix<T>& rhs);

    /**
     * @brief Element wise multiplication of a matrix by a vector.
     *
     * Returns a new matrix where each row is multiplied by rhs.
     */
    template <std::floating_point T>
    Matrix<T> hadamard(const Matrix<T>& lhs,
        const Vector<T>& rhs);

    /**
     * @brief In place element wise multiplication of a matrix by a vector.
     *
     * Each row of lhs is multiplied by rhs.
     */
    template<std::floating_point T>
    void hadamardInplace(Matrix<T>& lhs,
        const Vector<T>& rhs);

    /**
     * @brief Element wise square of a matrix.
     *
     * Returns a new matrix with each element squared.
     */
    template <std::floating_point T>
    Matrix<T> square(const Matrix<T>& m);

    /**
     * @brief In place element wise square of a matrix.
     */
    template <std::floating_point T>
    void squareInplace(Matrix<T>& m);

    /**
     * @brief Element wise square root of a matrix.
     *
     * Returns a new matrix with the square root applied to each element.
     */
    template <std::floating_point T>
    Matrix<T> sqrt(const Matrix<T>& m);

    /**
     * @brief In place element wise square root of a matrix.
     */
    template <std::floating_point T>
    void sqrtInplace(Matrix<T>& m);

    /**
     * @brief Convert a matrix to another value type.
     *
     * Returns a new matrix by static casting each element.
     */
    template <typename To, typename From>
    Matrix<To> convertMatrix(const Matrix<From>& A) noexcept;

} // namespace uv::core

#include "Functions.inl"