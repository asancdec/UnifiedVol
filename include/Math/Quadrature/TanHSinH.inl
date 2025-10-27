/**
* TanHSinH.inl
* Author: Alvaro Sanchez de Carlos
*/

#include <cmath>
#include <limits>

template<typename F>
double TanHSinH::integrateZeroToInf(F&& f) const noexcept
{
    // Machine Epsilon
	const long double eps{ std::numeric_limits<long double>::epsilon() };

	// Counter of how many terms fall below the relative error threshold
	unsigned int smallStreak{ 0U };

	// Running total of the right hand side compensated sum
	long double sumRight{ 0.0L };   

	// Compensation accumulator for lost low-order bits
	long double c{ 0.0L };     

	// Right side only
	for (const auto& node : nodes_)
	{
		// u(x) * w   (u includes the Jacobian; implement transformIntegrand(x,f))
		const long double term{ node.w * transformIntegrand(node.x, node.y, f) };

		// Guard against NaN tails
		if (!std::isfinite(term)) continue; 

		// Simple relative threshold vs current |sum|
		if (std::fabs(term) <= std::fabs(sumRight * eps))
		{	
			// Stop integrating if there are two insignificant consecutive evals
			if (++smallStreak >= 2) break;
		}
		else 
		{	
			// Reset the counter
			smallStreak = 0;
		}

		// Neumaier compensated summation
		const long double t{ sumRight + term };
		if (std::fabs(sumRight) >= std::fabs(term))
			c += (sumRight - t) + term;
		else
			c += (term - t) + sumRight;  
		sumRight = t;
	}

	// Add compensation term to the right-hand sum
	sumRight += c;

	// Running total of the left hand side compensated sum
	long double sumLeft{ 0.0L };

	// Reset variables
    smallStreak = 0U;
	c = 0.0L;

	// Left side only
	bool skipCenter{ true }; 
	for (const auto& node : nodes_)
	{
		if (skipCenter) 
		{	
			// Do not double count n=0
			skipCenter = false; 
			continue; 
		}

		// u(x) * w   (u includes the Jacobian; implement transformIntegrand(x,f))
		const long double term{ node.w * transformIntegrand(-node.x, 2.0L-node.y, f) };

		// Guard against NaN tails
		if (!std::isfinite(term)) continue;

		// Simple relative threshold vs current |sum|
		if (std::fabs(term) <= std::fabs(sumLeft * eps))
		{
			// Stop integrating if there are two insignificant consecutive evals
			if (++smallStreak >= 2) break;
		}
		else
		{
			// Reset the counter
			smallStreak = 0;
		}

		// Neumaier compensated summation
		const long double t{ sumLeft + term };
		if (std::fabs(sumLeft) >= std::fabs(term))
			c += (sumLeft - t) + term;   
		else
			c += (term - t) + sumLeft;
		sumLeft = t;
	}

	// Add compensation term to the left-hand sum
	sumLeft += c;

	return static_cast<double>(h_ * (sumLeft + sumRight));
}

template<typename F>
long double TanHSinH::transformIntegrand(long double x, long double y, F&& f) noexcept
{
	return 2.0L / (y * y) * static_cast<long double>(f((1.0L + x) / y));
}