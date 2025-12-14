// SPDX-License-Identifier: Apache-2.0
/*
 * File:        TanHSinH.inl
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

           
#include "Utils/IO/Log.hpp"     

#include <boost/math/special_functions/lambert_w.hpp>
#include <numbers>                                  
#include <cmath>                                  
#include <limits>                               
#include <iomanip>                                
#include <sstream>                               
#include <utility>
#include <iostream>
#include <bitset>

namespace uv::math
{
	template <std::size_t N>
	TanHSinH<N>::TanHSinH() :
		h_
		(
			// Define optimal step size using a heuristic rule
			boost::math::lambert_w0
			(
				Real(Real(2.0)) * std::numbers::pi_v<Real> * Real(N)
			)
			/ Real(N)
		)
	{	
		// Sanity checks at compile-time
		static_assert(N > 0,
			"TanHSinH<N>: N must be greater than zero.");
		static_assert(N % 2 == 0,
			"TanHSinH<N>: N must be even for unroll-by-2 integration.");

		// Calculate and store node values
		for (unsigned int n = 0; n < N; ++n)
		{
			nodes_[n] = generateNode(Real(n) * h_);
		}
	}

	template <std::size_t N>
	template<std::size_t M, typename F >
	std::array<Real, M> TanHSinH<N>::integrateZeroToInfMulti(F&& f) const noexcept
	{
		// Define numerical limit
		constexpr Real eps{ std::numeric_limits<Real>::epsilon() };

		// Accumulated sums
		std::array<Real, M>  sR0{ Real(Real(0.0))},
			sR1{ Real(Real(0.0))},
			sL0{ Real(Real(0.0))},
			sL1{ Real(Real(0.0))};

		// Bit sets to indicate when done
		std::bitset<M> actR0, actR1, actL0, actL1;
		actR0.set(); actR1.set(); actL0.set(); actL1.set();

		// Define forwarded callable
		auto&& func = std::forward<F>(f);

		// RHS: unrolled loop (process two nodes per iteration)
		for (std::size_t i = 0; i + 1 < N; i += 2)
		{
			// Check once per node whether any component is still active
			const bool anyA{ actR0.any() };
			const bool anyB{ actR1.any() };
			if (!anyA && !anyB) break;

			// Only evaluate func if at least one component needs it
			std::array<Real, M> ta;
			std::array<Real, M> tb;
			Real fa{ Real(Real(0.0)) }, fb{ Real(Real(0.0)) };
			if (anyA)
			{
				const Node& a{ nodes_[i] };
				ta = func(a.inputRight);
				fa = a.factorRight;
			};
			if (anyB)
			{
				const Node& b{ nodes_[i + 1] };
				tb = func(b.inputRight);
				fb = b.factorRight;
			};

			// Single pass over m
			for (std::size_t m = 0; m < M; ++m)
			{
				if (anyA && actR0.test(m))
				{
					// Early exit check (first)
					Real term{ fa * ta[m] };
					if (std::fabs(term) <= std::fabs(sR0[m] * eps)) actR0.reset(m);
					else sR0[m] += term;
				}
				if (anyB && actR1.test(m))
				{
					// Early exit check (second)
					Real term{ fb * tb[m] };
					if (std::fabs(term) <= std::fabs(sR1[m] * eps)) actR1.reset(m);
					else sR1[m] += term;
				}
			}
		 }

		// LHS: unrolled loop (process two nodes per iteration)
		for (std::size_t i = 1; i + 1 < N; i += 2)
		{
			// Check once per node whether any component is still active
			const bool anyA{ actL0.any() };
			const bool anyB{ actL1.any() };
			if (!anyA && !anyB) break;

			// Only evaluate func if at least one component needs it
			std::array<Real, M> ta;
			std::array<Real, M> tb;
			Real fa{ Real(Real(0.0)) }, fb{ Real(Real(0.0)) };
			if (anyA)
			{
				const Node& a{ nodes_[i] };
				ta = func(a.inputLeft);
				fa = a.factorLeft;
			};
			if (anyB)
			{
				const Node& b{ nodes_[i + 1] };
				tb = func(b.inputLeft);
				fb = b.factorLeft;
			};

			// Single pass over m
			for (std::size_t m = 0; m < M; ++m)
			{
				if (anyA && actL0.test(m))
				{
					// Early exit check (first)
					Real term{ fa * ta[m] };
					if (std::fabs(term) <= std::fabs(sL0[m] * eps)) actL0.reset(m);
					else sL0[m] += term;
				}
				if (anyB && actL1.test(m))
				{
					// Early exit check (second)
					Real term{ fb * tb[m] };
					if (std::fabs(term) <= std::fabs(sL1[m] * eps)) actL1.reset(m);
					else sL1[m] += term;
				}
			}
		}

		// Accumulate and return the total sum
		std::array<Real, M> out{};
		for (std::size_t m = 0; m < M; ++m) out[m] = sR0[m] + sR1[m] + sL0[m] + sL1[m];
		return out;
	}

	template <std::size_t N>
	template<typename F>
	Real TanHSinH<N>::integrateZeroToInf(F&& f) const noexcept
	{
		// Define numerical limit
		constexpr Real eps{ std::numeric_limits<Real>::epsilon() };

		// Accumulated sums
		Real sR0{ Real(Real(0.0)) }, sR1{ Real(Real(0.0)) }, sL0{ Real(Real(0.0)) }, sL1{ Real(Real(0.0)) };

		// Define forwarded callable
		auto&& func = std::forward<F>(f);

		// RHS: unrolled loop (process two nodes per iteration)
		for (std::size_t i = 0; i + 1 < N; i += 2)
		{	
			// n=0 is evaluated here
			const Node& a{ nodes_[i] };
			const Node& b{ nodes_[i + 1] };

			// u(x) * w  (first)
			const Real ta{ a.factorRight * func(a.inputRight) };

			// Early exit check (first)
			if (std::fabs(ta) <= std::fabs(sR0 * eps)) break;

			// Accumulate (first)
			sR0 += ta;

			// u(x) * w  (second)
			const Real tb{ b.factorRight * func(b.inputRight)};

			// Early exit check (second)
			if (std::fabs(tb) <= std::fabs(sR1 * eps)) break;

			// Accumulate (second)
			sR1 += tb;
		}

		// LHS: unrolled loop (process two nodes per iteration)
		for (std::size_t i = 1; i + 1 < N; i += 2)
		{
			const Node& a{ nodes_[i] };
			const Node& b{ nodes_[i + 1] };

			// u(x) * w  (first)
			const Real ta{ a.factorLeft * func(a.inputLeft) };

			// Early exit check (first)
			if (std::fabs(ta) <= std::fabs(sL0 * eps)) break;

			// Accumulate (first)
			sL0 += ta;

			// u(x) * w  (second)
			const Real tb{ b.factorLeft * func(b.inputLeft) };

			// Early exit check (second)
			if (std::fabs(tb) <= std::fabs(sL1 * eps)) break;

			// Accumulate (second)
			sL1 += tb;
		}

		// Accumulate and return the total sum
		return sR0 + sR1 + sL0 + sL1;
	}


	template <std::size_t N>
	void TanHSinH<N>::printGrid() const noexcept
	{
		constexpr int idxW = 6;
		constexpr int colW = 24;

		std::ostringstream oss;

		// Left align the title
		oss << std::left << "\nFixed Tanh-Sinh Grid\n";

		// Reset to default (right) for numeric columns
		oss << std::right;

		// Header
		oss << std::setw(colW) << "x_n (node)" << ' '
			<< std::setw(colW) << "w_n (weight)" << '\n';

		// Separator
		oss << std::string(idxW + 1 + colW + 1 + colW, '-') << '\n';

		// Body
		oss << std::scientific << std::setprecision(16);
		for (std::size_t i = 0; i < N; ++i)
		{
			oss << std::setw(idxW) << std::left << i
				<< std::setw(colW) << std::right << nodes_[i].x << ' '
				<< std::setw(colW) << std::right << nodes_[i].w << '\n';
		}

		UV_INFO(oss.str());
	}

	template <std::size_t N>
	TanHSinH<N>::Node TanHSinH<N>::generateNode(const Real nh) const noexcept
	{
		// Calculate q term
		const Real q{ std::exp(-std::numbers::pi_v<Real> * std::sinh(nh)) };

		// 1 / (1 + q)
		const Real qInv{ Real(Real(1.0)) / (Real(Real(1.0)) + q) };

		// Calculate y term
		const Real y{ Real(Real(2.0)) * q * qInv };

		// Calculate w
		const Real w{ qInv * y * std::numbers::pi_v<Real>*std::cosh(nh) };

		// Real(2.0) - y
		const Real twoMinusY{ Real(Real(2.0)) - y };

		// w * h
		const Real wh{w * h_};
		
		// Calculate and return Node struct
		return
		{
			w,                                                      // weight value
			y,                                                      // yn term
			(Real(Real(1.0)) - q) * qInv,                           // abscissas value
			wh * Real(Real(2.0)) / (y * y),                         // scaling term RHS
			twoMinusY / y,										    // transformed input RHS
			wh * Real(Real(2.0)) / (twoMinusY * twoMinusY),         // scaling term LHS
			y / twoMinusY 									        // transformed input LHS
		};
	}
}
