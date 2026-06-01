// SPDX-License-Identifier: Apache-2.0

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

    std::span<T> operator[](std::size_t) noexcept;

    std::span<const T> operator[](std::size_t) const noexcept;

    Matrix& operator+=(const Matrix&) noexcept;

    Matrix& operator-=(const Matrix&) noexcept;

    Matrix& operator+=(T) noexcept;

    Matrix& operator-=(T) noexcept;

    Matrix& operator*=(T) noexcept;

    Matrix& operator/=(T) noexcept;

    Matrix operator-() const noexcept;

    bool empty() const noexcept;

    std::size_t rows() const noexcept;

    std::size_t cols() const noexcept;

    template <std::floating_point U> Matrix<U> as() const noexcept;
};

template <std::floating_point T>
Matrix<T> operator+(Matrix<T>, const Matrix<T>&) noexcept;

template <std::floating_point T>
Matrix<T> operator-(Matrix<T>, const Matrix<T>&) noexcept;

template <std::floating_point T> Matrix<T> operator*(Matrix<T>, T) noexcept;

template <std::floating_point T> Matrix<T> operator/(Matrix<T>, T) noexcept;

template <std::floating_point T> Matrix<T> operator*(T, Matrix<T>) noexcept;

} // namespace uv::core

#include "Core/Detail/Matrix.inl"
