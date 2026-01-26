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
#include "Math/PDE/AHCache.hpp"
#include "Models/LocalVol/VarianceView.hpp"

#include <concepts>
#include <cstddef>
#include <span>
#include <array>

namespace uv::models::localvol
{
	template
	<
		std::floating_point T,
		std::size_t NT,
		std::size_t NX,
		class Interpolator
	>
	class Pricer
	{
	private:

		// ---------- Validate ----------

		// Size
		static_assert(NX >= 4, "NX must be >= 4");

		// Grids
		std::array<T, NX> xGrid_;
		std::array<T, NX> cInit_;
		std::array<T, NX> c_;

		// Views
		std::span<const T, NX - 2> xGridInner_;
		std::span<T, NX - 2> cInner_;

		// Cached variables and buffers
		math::pde::AHCache<T, NX> ahCache_{};

		// Interpolator
		Interpolator interpolator_{};

	public:

		// ---------- Typing ----------

		using interpolator_type = Interpolator;
		using derivatives_type = typename interpolator_type::derivatives_type;
		using evaluator_type   = typename interpolator_type::evaluator_type;


		Pricer() = delete;

		template <typename F>
		Pricer
		(
			F&& payoff,
			T xBound
		);

		Vector<T> priceNormalized
		(
			T tenor,
			std::span<const T> logKF,
			const VarianceView<T>& localVarView
		);

		Vector<T> price
		(
			T tenor,
			T forward,
			T discountFactor,
			std::span<const T> logKF,
			const VarianceView<T>& localVarView
		);

		Vector<T> price
		(
			Vector<T> tenor,
			T forward,
			T discountFactor,
			std::span<const T> logKF,
			const VarianceView<T>& localVarView
		);





	};

} // namespace uv::models::localvol

#include "Pricer.inl"