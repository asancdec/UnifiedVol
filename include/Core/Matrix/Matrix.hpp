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

#include "Core/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::core
{
	/**
	 * @brief Lightweight 2D matrix with contiguous row-major storage.
	 *
	 * The matrix stores all elements in a single contiguous buffer,
	 * enabling cache-friendly access and efficient SIMD auto-vectorization.
	 *
	 * All arithmetic operators are implemented as bounds-unchecked,
	 * in-place loops for maximum performance.
	 *
	 * @tparam T Floating-point element type.
	 */
	template <std::floating_point T>
	class Matrix
	{
	private:

		// ----------Member variables ----------

		std::size_t numRows_;
		std::size_t numCols_;
		Vector<T> data_;

	public:

		// ---------- Constructors ----------
		
		/** @brief Delete default constructor to avoid unwanted behaviour */
		Matrix() = delete;

		/**
		 * @brief Construct a matrix with given dimensions and initial value.
		 *
		 * Allocates a contiguous block of memory of size
		 * (numRows × numCols) and initializes all entries to the given value.
		 *
		 * @param numRows     Number of rows.
		 * @param numCols  Number of columns.
		 * @param val         Initial value for all elements (default = 0).
		 */
		Matrix(std::size_t numRows,
			std::size_t numCols,
			T val = T(0)) noexcept;

		// ---------- Unary operators ----------

		/** @brief Returns a mutable view of row @p i (no bounds checking). */
		std::span<T> operator[](std::size_t i) noexcept;

		/** @brief Returns a read-only view of row @p i (no bounds checking). */
		std::span<const T> operator[](std::size_t i) const noexcept;

		/** @brief In-place element-wise addition with another matrix. */
		Matrix& operator+=(const Matrix& rhs) noexcept;

		/** @brief In-place element-wise subtraction with another matrix. */
		Matrix& operator-=(const Matrix& rhs) noexcept;

		/** @brief In-place addition of a scalar to all elements. */
		Matrix& operator+=(T scalar) noexcept;

		/** @brief In-place subtraction of a scalar from all elements. */
		Matrix& operator-=(T scalar) noexcept;

		/** @brief In-place multiplication of all elements by a scalar. */
		Matrix& operator*=(T scalar) noexcept;

		/** @brief In-place division of all elements by a scalar. */
		Matrix& operator/=(T scalar) noexcept;

		/** @brief Returns a new matrix with all elements negated. */
		Matrix operator-() const noexcept;

		// ---------- Utilities ----------

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

		// ---------- Getters ----------

		/**
		 * @brief Number of rows in the matrix.
		 */
		std::size_t rows() const noexcept;

		/**
		 * @brief Number of columns in the matrix.
		 */
		std::size_t cols() const noexcept;

	};

	// ---------- Binary operators ----------

	/** @brief Returns the element-wise sum of two matrices. */
	template <std::floating_point T>
	Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs) noexcept;

	/** @brief Returns the element-wise difference of two matrices. */
	template <std::floating_point T>
	Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs) noexcept;

	/** @brief Returns a matrix scaled by a scalar (right multiplication). */
	template <std::floating_point T>
	Matrix<T> operator*(Matrix<T> lhs, T scalar) noexcept;

	/** @brief Returns a matrix divided by a scalar. */
	template <std::floating_point T>
	Matrix<T> operator/(Matrix<T> lhs, T scalar) noexcept;

	/** @brief Returns a matrix scaled by a scalar (left multiplication). */
	template <std::floating_point T>
	Matrix<T> operator*(T scalar, Matrix<T> rhs) noexcept;


} // namespace uv::core

#include "Matrix.inl"
