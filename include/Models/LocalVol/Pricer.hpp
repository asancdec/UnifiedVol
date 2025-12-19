// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.hpp
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

namespace uv::models::localvol
{
	template <std::floating_point T>
	class Pricer
	{ 
	private:


		Vector<T> r_;
		Vector<T> tGrid_;
		Vector<T> xGrid_;

		std::size_t nt_;
		std::size_t nx_;
		std::size_t nxMid_;

		T dt_;
		T dx_;

		Vector<T> initGuess_;
		Vector<T> xMidGrid_;

	public:

		Pricer() = delete;
		Pricer(std::span<const T> r,
			std::span<const T> tGrid,
			std::span<const T> xGrid);


		T price(T sigma) const;
	};
} // namespace uv::models::localvol

#include "Pricer.inl"