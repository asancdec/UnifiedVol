/*
* Types.hpp
* Author: Alvaro Sanchez de Carlos
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
     * @brief 2D dynamic container representing a numerical matrix.
     * 
     * @tparam T Element type.
     */
    template <typename T>
    using Matrix = std::vector<Vector<T>>;

    /**
     * @brief Complex number type with floating point real/imaginary components.
     *
     * @tparam T A floating point type (float, double, long double, etc.).
     */
    template<std::floating_point T>
    using Complex = std::complex<T>;
} // namespace uv