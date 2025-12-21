// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
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

namespace uv::math::integration
{
	template <std::floating_point T, std::size_t N, typename F>
	constexpr T trapezoidalWeighted(F&& f,
		const std::array<T, N>& y,
		const std::array<T, N>& x) noexcept
	{
		// ---------- Handle edge case ----------

		if constexpr (N < 2)
		{
			return T{ 0 };
		}
		else
		{
			// ---------- Integrate ----------

			// Prefer passing reference inside a loop
			auto&& fn = f;

			// First step
			T prevX{ x[0] };
			T prevFY{ static_cast<T>(fn(prevX)) * y[0] };
			T sum{ 0 };

			for (std::size_t i = 1; i < N; ++i)
			{
				// Extract
				const T curX{ x[i] };

				// Evaluate
				const T curFY{ static_cast<T>(fn(curX)) * y[i] };

				// Trapezoidal sum
				sum += (curX - prevX) * (prevFY + curFY) * T { 0.5 };

				// Update
				prevX = curX;
				prevFY = curFY;
			}

			return sum;
		}
	}

	template <std::floating_point T, std::size_t N>
	constexpr T trapezoidal(std::span< T, N> y,
		T dx) noexcept
	{
		// ---------- Handle edge case ----------

		if constexpr (N < 2)
		{
			return T{ 0 };
		}
		else
		{
			// ---------- Integrate ----------

			// Bounds
			T sum{ T{ 0.5 } * (y[0] + y[N - 1]) };

			// Inner domain
			for (std::size_t i = 1; i < N - 1; ++i)
			{
				sum += y[i];
			}

			return sum * dx;
		}  
	}

} // namespace uv::math::integration
