// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Matrix.hpp
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


namespace uv::core
{
    template <std::floating_point T>
    MatrixT<T>::MatrixT(const T val,
        const std::size_t numRows,
        const std::size_t numColumns) noexcept
        : numRows_(numRows),
        numColumns_(numColumns),
        data_(numRows * numColumns, val) {}

    template <std::floating_point T>
    std::size_t MatrixT<T>::rows() const noexcept
    {
        return numRows_;
    }

    template <std::floating_point T>

    std::span<T> MatrixT<T>::operator[](std::size_t i) noexcept
    {      
        return { data_.data() + i * numColumns_, numColumns_ };
    }

    template <std::floating_point T>
    std::span<const T> MatrixT<T>::operator[](std::size_t i) const noexcept 
    {
        return { data_.data() + i * numColumns_, numColumns_ };
    }



    //template <std::floating_point T>
    //void MatrixT<T>::print() const noexcept
    //{
    //    uv::utils::printMatrix(
    //        /*title=*/title,
    //        /*header=*/header,
    //        /*rowLabels=*/tenors_,
    //        /*M=*/volMatrix_,
    //        /*headerPrec=*/headerPrec,
    //        /*rowLabelPrec=*/2,
    //        /*valuePrec=*/valuePrec
    //    );
    //}


} // namespace uv::core