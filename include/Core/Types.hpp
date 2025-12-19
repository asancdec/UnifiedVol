// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Types.hpp
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


#pragma once

#include <complex>
#include <type_traits>
#include <vector>
#include <concepts>

namespace uv
{
    /**
     * @brief Real number type used throughout the library.
     *
     * CHANGE HERE if you want global numeric precision:
     *     - `double` for speed
     *     - `long double` for precision
     */
    using Real = double;
    static_assert(std::is_floating_point_v<Real>,
        "uv::Real must be a floating point type");

    /**
     * @brief 1D dynamic contiguous container (vector) of values.
     *
     * @tparam T Element type.
     */
    template <typename T>
    using Vector = std::vector<T>;

    /**
     * @brief Complex number type with floating point real/imaginary components.
     *
     * @tparam T A floating point type (float, double, long double, etc.).
     */
    template<std::floating_point T>
    using Complex = std::complex<T>;
} // namespace uv