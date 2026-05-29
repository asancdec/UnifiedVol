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

#include "Base/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::core
{

template <std::floating_point T> class Matrix
{
  private:
    std::size_t numRows_;
    std::size_t numCols_;
    Vector<T> data_;

  public:
    Matrix() = delete;

    Matrix(std::size_t numRows, std::size_t numCols, T val = 0.0) noexcept;

    std::span<T> operator[](std::size_t i) noexcept;

    std::span<const T> operator[](std::size_t i) const noexcept;

    Matrix& operator+=(const Matrix& rhs) noexcept;

    Matrix& operator-=(const Matrix& rhs) noexcept;

    Matrix& operator+=(T scalar) noexcept;

    Matrix& operator-=(T scalar) noexcept;

    Matrix& operator*=(T scalar) noexcept;

    Matrix& operator/=(T scalar) noexcept;

    Matrix operator-() const noexcept;

    bool empty() const noexcept;

    std::size_t rows() const noexcept;

    std::size_t cols() const noexcept;

    template <std::floating_point U> Matrix<U> as() const noexcept;
};

template <std::floating_point T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs) noexcept;

template <std::floating_point T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs) noexcept;

template <std::floating_point T> Matrix<T> operator*(Matrix<T> lhs, T scalar) noexcept;

template <std::floating_point T> Matrix<T> operator/(Matrix<T> lhs, T scalar) noexcept;

template <std::floating_point T> Matrix<T> operator*(T scalar, Matrix<T> rhs) noexcept;

} // namespace uv::core

#include "Core/Detail/Matrix.inl"
