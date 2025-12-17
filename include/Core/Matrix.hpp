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

#pragma once

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::core
{
	/**
	 * @brief Lightweight 2D matrix with contiguous memory storage.
	 * Internally, data is stored in a single contiguous 1D buffer 
	 * in row-major order, enabling:
	 *
	 * - Cache-friendly access patterns
	 * - Efficient SIMD auto-vectorization
	 * - Zero-copy row views via std::span
	 *
	 * The class provides bounds-unchecked row access through operator[],
	 * returning a std::span over the requested row.
	 *
	 * @tparam T Floating-point value type.
	 */
	template <std::floating_point T>
	class Matrix
	{
	private:

		std::size_t numRows_{ 0 };
		std::size_t numColumns_{ 0 };
		Vector<T> data_{};

	public:

		Matrix() = delete;

		/**
		 * @brief Construct a matrix with given dimensions and initial value.
		 *
		 * Allocates a contiguous block of memory of size
		 * (numRows × numColumns) and initializes all entries to the given value.
		 *
		 * @param numRows     Number of rows.
		 * @param numColumns  Number of columns.
		 * @param val         Initial value for all elements (default = 0).
		 */
		Matrix(std::size_t numRows,
		std::size_t numColumns,
		T val = T(0)) noexcept;

		/**
		 * @brief Access a matrix row as a contiguous view.
		 *
		 * Returns a std::span representing the i-th row of the matrix.
		 * The returned span provides direct access to the underlying data
		 * without copying.
		 *
		 * @param i Row index (0-based).
		 * @return Mutable span over the row elements.
		 *
		 * @warning No bounds checking is performed.
		 */
		std::span<T> operator[](std::size_t i) noexcept;

		/**
		 * @brief Access a matrix row as a read-only contiguous view.
		 *
		 * Returns a const std::span representing the i-th row of the matrix.
		 *
		 * @param i Row index (0-based).
		 * @return Read-only span over the row elements.
		 *
		 * @warning No bounds checking is performed.
		 */
		std::span<const T> operator[](std::size_t i) const noexcept;

		/**
		 * @brief Print the matrix to standard output.
		 *
		 * Prints the matrix with optional formatting precision for values.
		 * Intended for debugging and diagnostics only.
		 *
		 * @param valuePrec Number of decimal places for matrix values.
		 */
		void print(unsigned int valuePrec = 4) const noexcept;

		/**
		 * @brief Check whether the matrix is empty.
		 *
		 * A matrix is considered empty if it has zero rows or zero columns.
		 *
		 * @return True if the matrix is empty, false otherwise.
		 */
		bool empty() const noexcept;

		/**
		 * @brief Number of rows in the matrix.
		 */
		std::size_t rows() const noexcept;

		/**
		 * @brief Number of columns in the matrix.
		 */
		std::size_t cols() const noexcept;

	};

} // namespace uv::core

#include "Matrix.inl"
