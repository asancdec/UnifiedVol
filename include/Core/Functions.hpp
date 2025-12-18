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

#include <span>
#include <cstddef>

namespace uv::core
{
    /**
     * @brief Generate an evenly spaced 1D grid.
     *
     * Returns a vector of values from bound1 to bound2 (inclusive).
     */
    Vector<Real> generateGrid(Real bound1,
        Real bound2,
        size_t steps
    ) noexcept;

    /**
     * @brief Generate a sequence of consecutive values.
     *
     * Returns a vector of length n starting from start.
     */
    template <typename T>
    Vector<T> makeSequence(std::size_t n, T start) noexcept;

    /**
     * @brief Compute forward differences of a vector.
     *
     * Returns v[i+1] - v[i] for all valid indices.
     */
    Vector<Real> diff(const Vector<Real>& v);

    /**
     * @brief Multiply all elements of a vector by a scalar.
     */
    Vector<Real> multiply(const Vector<Real>& v,
        Real x) noexcept;

    /**
     * @brief Compute element wise reciprocal of a vector.
     */
    Vector<Real> reciprocal(const Vector<Real>& v) noexcept;

    /**
     * @brief Element wise multiplication of two vectors.
     */
    Vector<Real> hadamard(std::span<const Real> a,
        std::span<const Real> b);

    /**
     * @brief Return the minimum element of a vector.
     *
     * The input must be non empty.
     */
    template<typename T>
    T minValue(std::span<const T> x);

    /**
     * @brief Return the maximum element of a vector.
     *
     * The input must be non empty.
     */
    template<typename T>
    T maxValue(std::span<const T> x);

    /**
     * @brief Convert a vector to another value type.
     *
     * Returns a new vector with each element static casted.
     */
    template <typename To, typename From>
    Vector<To> convertVector(const Vector<From>& x) noexcept;

} // namespace uv::core

#include "Functions.inl"