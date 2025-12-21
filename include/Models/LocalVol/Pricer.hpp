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

#include "Core/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>
#include <array>

namespace uv::models::localvol
{
	template 
	<
		std::floating_point T, 
		std::size_t nT, 
		std::size_t nX,
		typename Fn
	>
	class Pricer
	{ 
	private:

		// Chang-Cooper coefficients
		std::array<T, nX - 1 > B_{};
		std::array<T, nX - 1 > C_{};

		// Payoff data
		T t_;
		T r_;
		T F_;
		T K_;
		Fn payoff_;

		// Grid data
		std::array<T, nX> xGrid_;
		T dt_;
		T dx_;
		std::array<T, nX> pdfGrid_;

		// Private functions
		T normalizeForward_(T F) const noexcept;
		
	public:

		Pricer() = delete;
		Pricer(T t,
			T r,
			T F,
			T K,
			Fn payoff,
			const std::array<T, nX>& xGrid);

		T price(T localVar);
	};
} // namespace uv::models::localvol

#include "Pricer.inl"