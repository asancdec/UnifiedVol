// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Matrix.inl
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


#include "Utils/IO/Functions.hpp"
#include "Core/Functions.hpp"

#include <numeric>
#include <string_view>

namespace uv::core
{
    template <std::floating_point T>
    Matrix<T>::Matrix(const std::size_t numRows,
        std::size_t numColumns,
        T val) noexcept
        : numRows_(numRows),
        numColumns_(numColumns),
        data_(numRows * numColumns, val) {}

    template <std::floating_point T>
    std::span<T> Matrix<T>::operator[](std::size_t i) noexcept
    {      
        return { data_.data() + i * numColumns_, numColumns_ };
    }

    template <std::floating_point T>
    std::span<const T> Matrix<T>::operator[](std::size_t i) const noexcept 
    {
        return { data_.data() + i * numColumns_, numColumns_ };
    }

    template <std::floating_point T>
    void Matrix<T>::print(unsigned int valuePrec) const noexcept
    {
        utils::printMatrix(
            /*title=*/"",
            /*header=*/makeSequence<T>(numColumns_, 1),
            /*rowLabels=*/makeSequence<T>(numRows_, 1),
            /*M=*/*this,
            /*headerPrec=*/1,
            /*rowLabelPrec=*/1,
            /*valuePrec=*/valuePrec
        );
    }

    template <std::floating_point T>
    bool Matrix<T>::empty() const noexcept
    {
        return numRows_ == 0 || numColumns_ == 0;
    }

    template <std::floating_point T>
    std::size_t Matrix<T>::rows() const noexcept
    {
        return numRows_;
    }

    template <std::floating_point T>
    std::size_t Matrix<T>::cols() const noexcept
    {
        return numColumns_;
    }

} // namespace uv::core