/**
* TanHSinH.inl
* Author: Alvaro Sanchez de Carlos
*/

#include <cmath>
#include <limits>

template<typename F>
long double TanHSinH::integrateZeroToInf(F&& f) const noexcept
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
		if (!std::isfinite(term)) continue;

		// Check if calculated term is below relative threshold
		if (std::fabs(term) <= std::fabs(sumRight * std::numeric_limits<long double>::epsilon()))
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
	for (std::size_t i = 1; i < nodes_.size(); ++i)
	{
		// Define node
		const Node& node{ nodes_[i] };

		// u(x) * w 
		const long double term{ node.w * transformIntegrand(-node.x, 2.0L - node.y, f) };

		// Guard against NaN
		if (!std::isfinite(term)) continue;

		// Check if calculated term is below relative threshold
		if (std::fabs(term) <= std::fabs(sumLeft * std::numeric_limits<long double>::epsilon()))
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

template<typename F>
long double TanHSinH::transformIntegrand(long double x, long double y, F&& f) noexcept
{
	return 2.0L / (y * y) * static_cast<long double>(f((1.0L + x) / y));
}