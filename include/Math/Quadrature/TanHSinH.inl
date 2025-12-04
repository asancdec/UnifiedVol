/**
* TanHSinH.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Quadrature/TanHSinH.hpp"              
#include "Utils/Log.hpp"                             

#include <boost/math/special_functions/lambert_w.hpp>
#include <numbers>                                  
#include <cmath>                                  
#include <limits>                               
#include <iomanip>                                
#include <sstream>                               
#include <utility>
#include <iostream>
#include <bitset>

namespace uv
{
	template <::std::size_t N>
	TanHSinH<N>::TanHSinH() :
		h_
		(
			// Define optimal step size using a heuristic rule
			boost::math::lambert_w0
			(
				2.0L * ::std::numbers::pi_v<long double> *static_cast<long double>(N)
			)
			/ static_cast<long double>(N)
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
			nodes_[n] = generateNode(static_cast<long double>(n) * h_);
		}
	}

	template <::std::size_t N>
	template<::std::size_t M, typename F >
	::std::array<long double, M> TanHSinH<N>::integrateZeroToInfMulti(F&& f) const noexcept
	{
		// Define numerical limit
		constexpr long double eps{ ::std::numeric_limits<long double>::epsilon() };

		// Accumulated sums
		::std::array<long double, M>  sR0{ 0.0L }, sR1{ 0.0L }, sL0{ 0.0L }, sL1{ 0.0L };

		// Bit sets to indicate when done
		::std::bitset<M> actR0, actR1, actL0, actL1;
		actR0.set(); actR1.set(); actL0.set(); actL1.set();

		// Define forwarded callable
		auto&& func = std::forward<F>(f);

		// RHS: unrolled loop (process two nodes per iteration)
		for (::std::size_t i = 0; i + 1 < N; i += 2)
		{
			// Check once per node whether any component is still active
			const bool anyA{ actR0.any() };
			const bool anyB{ actR1.any() };
			if (!anyA && !anyB) break;

			// Only evaluate func if at least one component needs it
			::std::array<long double, M> ta;
			::std::array<long double, M> tb;
			long double fa{ 0.0L }, fb{ 0.0L };
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
			for (::std::size_t m = 0; m < M; ++m)
			{
				if (anyA && actR0.test(m))
				{
					// Early exit check (first)
					long double term{ fa * ta[m] };
					if (::std::fabs(term) <= ::std::fabs(sR0[m] * eps)) actR0.reset(m);
					else sR0[m] += term;
				}
				if (anyB && actR1.test(m))
				{
					// Early exit check (second)
					long double term{ fb * tb[m] };
					if (::std::fabs(term) <= ::std::fabs(sR1[m] * eps)) actR1.reset(m);
					else sR1[m] += term;
				}
			}
		 }

		// LHS: unrolled loop (process two nodes per iteration)
		for (::std::size_t i = 1; i + 1 < N; i += 2)
		{
			// Check once per node whether any component is still active
			const bool anyA{ actL0.any() };
			const bool anyB{ actL1.any() };
			if (!anyA && !anyB) break;

			// Only evaluate func if at least one component needs it
			::std::array<long double, M> ta;
			::std::array<long double, M> tb;
			long double fa{0.0L}, fb{ 0.0L };
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
			for (::std::size_t m = 0; m < M; ++m)
			{
				if (anyA && actL0.test(m))
				{
					// Early exit check (first)
					long double term{ fa * ta[m] };
					if (::std::fabs(term) <= ::std::fabs(sL0[m] * eps)) actL0.reset(m);
					else sL0[m] += term;
				}
				if (anyB && actL1.test(m))
				{
					// Early exit check (second)
					long double term{ fb * tb[m] };
					if (::std::fabs(term) <= ::std::fabs(sL1[m] * eps)) actL1.reset(m);
					else sL1[m] += term;
				}
			}
		}

		// Accumulate and return the total sum
		::std::array<long double, M> out{};
		for (std::size_t m = 0; m < M; ++m) out[m] = sR0[m] + sR1[m] + sL0[m] + sL1[m];
		return out;
	}

	template <::std::size_t N>
	template<typename F>
	long double TanHSinH<N>::integrateZeroToInf(F&& f) const noexcept
	{
		// Define numerical limit
		constexpr long double eps{ ::std::numeric_limits<long double>::epsilon() };

		// Accumulated sums
		long double sR0{ 0.0L }, sR1{ 0.0L }, sL0{ 0.0L }, sL1{ 0.0L };

		// Define forwarded callable
		auto&& func = std::forward<F>(f);

		// RHS: unrolled loop (process two nodes per iteration)
		for (::std::size_t i = 0; i + 1 < N; i += 2)
		{	
			// n=0 is evaluated here
			const Node& a{ nodes_[i] };
			const Node& b{ nodes_[i + 1] };

			// u(x) * w  (first)
			const long double ta{ a.factorRight * static_cast<long double>(func(a.inputRight)) };

			// Early exit check (first)
			if (::std::fabs(ta) <= ::std::fabs(sR0 * eps)) break;

			// Accumulate (first)
			sR0 += ta;

			// u(x) * w  (second)
			const long double tb{ b.factorRight * static_cast<long double>(func(b.inputRight)) };

			// Early exit check (second)
			if (::std::fabs(tb) <= ::std::fabs(sR1 * eps)) break;

			// Accumulate (second)
			sR1 += tb;
		}

		// LHS: unrolled loop (process two nodes per iteration)
		for (::std::size_t i = 1; i + 1 < N; i += 2)
		{
			const Node& a{ nodes_[i] };
			const Node& b{ nodes_[i + 1] };

			// u(x) * w  (first)
			const long double ta{ a.factorLeft * static_cast<long double>(func(a.inputLeft)) };

			// Early exit check (first)
			if (::std::fabs(ta) <= ::std::fabs(sL0 * eps)) break;

			// Accumulate (first)
			sL0 += ta;

			// u(x) * w  (second)
			const long double tb{ b.factorLeft * static_cast<long double>(func(b.inputLeft)) };

			// Early exit check (second)
			if (::std::fabs(tb) <= ::std::fabs(sL1 * eps)) break;

			// Accumulate (second)
			sL1 += tb;
		}

		// Accumulate and return the total sum
		return sR0 + sR1 + sL0 + sL1;
	}


	template <::std::size_t N>
	void TanHSinH<N>::printGrid() const noexcept
	{
		constexpr int idxW = 6;
		constexpr int colW = 24;

		::std::ostringstream oss;

		// Left align the title
		oss << ::std::left << "\nFixed Tanh-Sinh Grid\n";

		// Reset to default (right) for numeric columns
		oss << ::std::right;

		// Header
		oss << ::std::setw(colW) << "x_n (node)" << ' '
			<< ::std::setw(colW) << "w_n (weight)" << '\n';

		// Separator
		oss << ::std::string(idxW + 1 + colW + 1 + colW, '-') << '\n';

		// Body
		oss << ::std::scientific << ::std::setprecision(16);
		for (::std::size_t i = 0; i < N; ++i)
		{
			oss << ::std::setw(idxW) << ::std::left << i
				<< ::std::setw(colW) << ::std::right << nodes_[i].x << ' '
				<< ::std::setw(colW) << ::std::right << nodes_[i].w << '\n';
		}

		UV_INFO(oss.str());
	}

	template <::std::size_t N>
	TanHSinH<N>::Node TanHSinH<N>::generateNode(const long double nh) const noexcept
	{
		// Calculate q term
		const long double q{ ::std::exp(-::std::numbers::pi_v<long double> * ::std::sinh(nh)) };

		// 1 / (1 + q)
		const long double qInv{ 1.0L / (1.0L + q) };

		// Calculate y term
		const long double y{ 2.0L * q * qInv };

		// Calculate w
		const long double w{ qInv * y * ::std::numbers::pi_v<long double>*::std::cosh(nh) };

		// 2.0 - y
		const long double twoMinusY{ 2.0L - y };

		// w * h
		const long double wh{w * h_};
		
		// Calculate and return Node struct
		return
		{
			w,                                                                // weight value
			y,                                                                // yn term
			(1.0L - q) * qInv,                                                // abscissas value
			wh * 2.0L / (y * y),                                              // scaling term RHS
			twoMinusY / y,												      // transformed input RHS
			wh * 2.0L / (twoMinusY * twoMinusY),                              // scaling term LHS
			y / twoMinusY 												      // transformed input LHS
		};
	}
}
