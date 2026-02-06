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

#include <numeric>
#include <string_view>

namespace uv::core
{
template <std::floating_point T>
Matrix<T>::Matrix(const std::size_t numRows, std::size_t numCols, T val) noexcept
    : numRows_(numRows),
      numCols_(numCols),
      data_(numRows * numCols, val)
{
}

template <std::floating_point T>
std::span<T> Matrix<T>::operator[](std::size_t i) noexcept
{
    return {data_.data() + i * numCols_, numCols_};
}

template <std::floating_point T>
std::span<const T> Matrix<T>::operator[](std::size_t i) const noexcept
{
    return {data_.data() + i * numCols_, numCols_};
}

template <std::floating_point T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) noexcept
{
    T* ptrThis{data_.data()};
    const T* ptrRhs{rhs.data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptrThis[i] += ptrRhs[i];
    }

    return *this;
}

template <std::floating_point T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) noexcept
{
    T* ptrThis{data_.data()};
    const T* ptrRhs{rhs.data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptrThis[i] -= ptrRhs[i];
    }

    return *this;
}

template <std::floating_point T> Matrix<T>& Matrix<T>::operator+=(T scalar) noexcept
{
    T* ptr{data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptr[i] += scalar;
    }

    return *this;
}

template <std::floating_point T> Matrix<T>& Matrix<T>::operator-=(T scalar) noexcept
{
    T* ptr{data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptr[i] -= scalar;
    }

    return *this;
}

template <std::floating_point T> Matrix<T>& Matrix<T>::operator*=(T scalar) noexcept
{
    T* ptr{data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptr[i] *= scalar;
    }

    return *this;
}

template <std::floating_point T> Matrix<T>& Matrix<T>::operator/=(T scalar) noexcept
{
    T* ptr{data_.data()};

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        ptr[i] /= scalar;
    }

    return *this;
}

template <std::floating_point T> Matrix<T> Matrix<T>::operator-() const noexcept
{
    Matrix<T> result(*this);

    T* ptr{result.data_.data()};

    for (std::size_t i = 0; i < result.data_.size(); ++i)
    {
        ptr[i] = -ptr[i];
    }

    return result;
}

template <std::floating_point T> bool Matrix<T>::empty() const noexcept
{
    return numRows_ == 0 || numCols_ == 0;
}

template <std::floating_point T> std::size_t Matrix<T>::rows() const noexcept
{
    return numRows_;
}

template <std::floating_point T> std::size_t Matrix<T>::cols() const noexcept
{
    return numCols_;
}

template <std::floating_point T>
template <std::floating_point U>
Matrix<U> Matrix<T>::as() const noexcept
{
    Matrix<U> out(numRows_, numCols_);

    for (std::size_t i{0}; i < numRows_; ++i)
    {
        std::span<const T> rowFrom{(*this)[i]};
        std::span<U> rowTo{out[i]};

        for (std::size_t j{0}; j < numCols_; ++j)
        {
            rowTo[j] = static_cast<U>(rowFrom[j]);
        }
    }

    return out;
}

template <std::floating_point T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs) noexcept
{
    lhs += rhs;
    return lhs;
}

template <std::floating_point T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs) noexcept
{
    lhs -= rhs;
    return lhs;
}

template <std::floating_point T> Matrix<T> operator*(Matrix<T> lhs, T scalar) noexcept
{
    lhs *= scalar;
    return lhs;
}

template <std::floating_point T> Matrix<T> operator/(Matrix<T> lhs, T scalar) noexcept
{
    lhs /= scalar;
    return lhs;
}

template <std::floating_point T> Matrix<T> operator*(T scalar, Matrix<T> rhs) noexcept
{
    rhs *= scalar;
    return rhs;
}

} // namespace uv::core