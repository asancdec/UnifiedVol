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

#include <algorithm>
#include <numeric>
#include <span>

namespace uv::core
{
    template <typename T>
    Vector<T> makeSequence(std::size_t n, 
        T start) noexcept
    {
        Vector<T> v(n);
        std::iota(v.begin(), v.end(), start);
        return v;
    }

    template<typename T>
    T minValue(std::span<const T> x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "minValue: input vector is empty"
        );

        return *std::min_element(x.begin(), x.end());
    }

    template<typename T>
    T maxValue(std::span<const T> x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "maxValue: input vector is empty"
        );

        return *std::max_element(x.begin(), x.end());
    }


    template <typename To, typename From>
    Vector<To> convertVector(const Vector<From>& x) noexcept
    {
        Vector<To> out;
        out.reserve(x.size());

        for (const From& v : x)
            out.push_back(static_cast<To>(v));

        return out;
    }

    template <typename To, typename From>
    Matrix<To> convertMatrix(const Matrix<From>& A) noexcept
    {
        std::size_t nRows{ A.rows() };
        std::size_t nCols{ A.cols() };

        core::Matrix<To> out(nRows, nCols);
 
        for (std::size_t i = 0; i < nRows; ++i)
        {
            std::span<const From> inRow{ A[i] };
            std::span<To> outRow{ out[i] };

            for (std::size_t j = 0; j < nCols; ++j)
            {
                outRow[j] = static_cast<To>(inRow[j]);
            }
        }
        return out;
    }

} // namespace uv::core