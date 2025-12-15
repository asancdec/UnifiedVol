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

#include "Utils/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::core
{
	template <std::floating_point T>
	class MatrixT
	{
	private:

		std::size_t numRows_{ 0 };
		std::size_t numColumns_{ 0 };
		Vector<T> data_{};

	public:

		MatrixT(T val,
		std::size_t numRows,
		std::size_t numColumns) noexcept;


		std::size_t rows() const noexcept;


		std::span<T> operator[](std::size_t i) noexcept;
		std::span<const T> operator[](std::size_t i) const noexcept;

		void print() const noexcept;

	};

} // namespace uv::core

#include "Matrix.inl"