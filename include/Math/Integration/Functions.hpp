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

#include <array>
#include <concepts>
#include <cstddef>

namespace uv::math::integration
{
	/**
	 * @brief Trapezoidal integration of a weighted function f(x) * y(x) on a non-uniform grid.
	 *
	 * Computes the integral of the pointwise product f(x[i]) * y[i] with respect to x
	 * using the trapezoidal rule on the provided grid `x` (assumed ordered).
	 *
	 * @tparam T Floating-point type.
	 * @tparam N Number of grid points.
	 * @tparam F Callable type, invoked as f(x).
	 * @param f Weighting function evaluated at x.
	 * @param y Samples (e.g., pdf values) on the grid.
	 * @param x Grid points (not necessarily equally spaced).
	 * @return Approximate integral of f(x) * y(x) dx (returns 0 if N < 2).
	 */
	template <std::floating_point T, std::size_t N, typename F>
	constexpr T trapezoidalWeighted(F&& f,
		const std::array<T, N>& y,
		const std::array<T, N>& x) noexcept;

	/**
	 * @brief Trapezoidal integration on a uniform grid and with a pre-calculated
	 * grid of y values
	 *
	 * Computes the integral of sampled values `y` assuming constant spacing `dx`:
	 * `dx * (0.5*y[0] + y[1] + ... + y[N-2] + 0.5*y[N-1])`.
	 *
	 * @param y  Samples on an evenly spaced grid (size N).
	 * @param dx Constant grid spacing.
	 * @return Approximate integral (returns 0 if N < 2).
	 */
	template <std::floating_point T, std::size_t N>
	constexpr T trapezoidal(std::span<T, N> y,
		T dx) noexcept;

} // namespace uv::math::integration

#include "Functions.inl"
