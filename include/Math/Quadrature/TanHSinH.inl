/**
* TanHSinH.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Quadrature/TanHSinH.hpp"              
#include "Errors/Errors.hpp"                       
#include "Utils/Log.hpp"                             

#include <boost/math/special_functions/lambert_w.hpp>
#include <numbers>                                  
#include <cmath>                                  
#include <limits>                               
#include <iomanip>                                
#include <sstream>                               

namespace uv
{
	template <::std::size_t N>
	TanHSinH<N>::TanHSinH() :
		h_
		(
			// Define optimal step size using a heuristic rule
			boost::math::lambert_w0
			(
				2.0L * ::std::numbers::pi_v<long double> * static_cast<long double>(N)
			)
			/ static_cast<long double>(N)
		)
	{
		// Calculate and store node values
		for (unsigned int n = 0; n < N; ++n)
		{
			nodes_[n] = generateNode(static_cast<long double>(n) * h_);
		}
	}

	template <::std::size_t N>
	TanHSinH<N>::Node TanHSinH<N>::generateNode(const long double nh) noexcept
	{
		// Calculate qn term
		long double qn{ ::std::exp(-::std::numbers::pi_v<long double> * ::std::sinh(nh)) };

		// Precompute repetitive calculations
		long double qnInv{ 1.0L / (1.0L + qn) };

		// Calculate yn term
		long double yn{ (2.0L * qn * qnInv) };

		// Calculate and return Node struct
		return
		{
			yn,                                                                 // yn term
			1.0L - yn,                                                          // Abscissas value
			qnInv * yn * ::std::numbers::pi_v<long double> *::std::cosh(nh)        // Weight value
		};
	}

	template <::std::size_t N>
	void TanHSinH<N>::printGrid() const noexcept
	{
		constexpr int idx_w = 6;
		constexpr int col_w = 24;

		::std::ostringstream oss;

		// Left align the title
		oss << ::std::left << "\nFixed Tanh-Sinh Grid\n";

		// Reset to default (right) for numeric columns
		oss << ::std::right;

		// Header
		oss << ::std::setw(col_w) << "x_n (node)" << ' '
			<< ::std::setw(col_w) << "w_n (weight)" << '\n';

		// Separator
		oss << ::std::string(idx_w + 1 + col_w + 1 + col_w, '-') << '\n';

		// Body
		oss << ::std::scientific << ::std::setprecision(16);
		for (::std::size_t n = 0; n < N; ++n)
		{
			oss << ::std::setw(idx_w) << ::std::left << n
				<< ::std::setw(col_w) << ::std::right << nodes_[n].x << ' '
				<< ::std::setw(col_w) << ::std::right << nodes_[n].w << '\n';
		}

		UV_INFO(oss.str());
	}

	template <::std::size_t N>
	template<typename F>
	long double TanHSinH<N>::integrateZeroToInf(F&& f) const noexcept
	{
		// Counter of the number of terms that fall below evaluation threshold
		unsigned int smallStreak{ 0U };

		// Accumulated sum of the right hand area
		long double sumRight{ 0.0L };

		// Right-hand side
		for (const auto& node : nodes_)
		{
			// u(x) * w 
			const long double term{ node.w * transformIntegrand(node.x, node.y, f) };

			// Guard against NaN
			if (!::std::isfinite(term)) continue;

			// Check if calculated term is below relative threshold
			if (::std::fabs(term) <= ::std::fabs(sumRight * ::std::numeric_limits<long double>::epsilon()))
			{
				// Stop integrating if two consecutive streaks
				if (++smallStreak >= 2U)
				{
					smallStreak = 0U;
					break;
				}
			}
			else
			{
				smallStreak = 0;
			}

			// Accumulate to current sum
			sumRight += term;
		}

		// Accumulated sum of the left hand area
		long double sumLeft{ 0.0L };

		// Left-hand side (do not double count n=0)
		for (::std::size_t i = 1; i < N; ++i)
		{
			// Define node
			const Node& node{ nodes_[i] };

			// u(x) * w 
			const long double term{ node.w * TanHSinH<N>::transformIntegrand(-node.x, 2.0L - node.y, f) };

			// Guard against NaN
			if (!::std::isfinite(term)) continue;

			// Check if calculated term is below relative threshold
			if (::std::fabs(term) <= ::std::fabs(sumLeft * ::std::numeric_limits<long double>::epsilon()))
			{
				// Stop integrating if two consecutive streaks
				if (++smallStreak >= 2U) break;
			}
			else
			{
				smallStreak = 0;
			}

			// Accumulate to current sum
			sumLeft += term;
		}

		return h_ * (sumLeft + sumRight);
	}

	template <::std::size_t N>
	template<typename F>
	long double TanHSinH<N>::transformIntegrand(long double x, long double y, F&& f) noexcept
	{
		return 2.0L / (y * y) * static_cast<long double>(f((1.0L + x) / y));
	}
}
