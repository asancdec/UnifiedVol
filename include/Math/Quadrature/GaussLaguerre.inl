/**
* GaussLaguerre.inl
* Author: Alvaro Sanchez de Carlos
*/

#include <cmath>

template<typename F>
double GaussLaguerre::eval(F&& f) const
{

	long double sum{ 0.0L };   // Running total of the compensated sum
	long double c{ 0.0L };     // Compensation accumulator for lost low-order bits

	for (int i = 0; i < N_; ++i)
	{	
		// Fused multiply-add: computes sum + wk[i]*f(xk[i]) with a single rounding
		const long double t{ std::fma(wk_[i], static_cast<long double>(f(static_cast<double>(xk_[i]))), sum) };

		// Increment actually added to sum
		// Any bits lost in this step are accounted for in the c-update below.
		const long double term{ t - sum };

		// Neumaier update
		if (std::fabs(sum) >= std::fabs(term))
			c += (sum - t) + term;   // lost bits when |sum| ≥ |term|
		else
			c += (term - t) + sum;   // lost bits when |term| > |sum|

		sum = t; // advance running total
	}
	return static_cast<double>(sum + c);
}

template<typename F>
double GaussLaguerre::evalUnweighted(F&& f) const
{
	long double sum{ 0.0L };   // Running total of the compensated sum
	long double c{ 0.0L };     // Compensation accumulator for lost low-order bits

	for (int i = 0; i < N_; ++i)
	{
		// Unweighting transform:
		//   ∫₀^∞ f(x) dx = ∫₀^∞ x^α e^{-x} [ e^{x} x^{-α} f(x) ] dx
		// so per-node term is e^{+x_i} * x_i^{−α} * f(x_i)
		const long double term{ std::exp(xk_[i]) * std::pow(xk_[i], -alpha_) 
			* static_cast<long double>(f(static_cast<double>(xk_[i])))};

		// Fused multiply-add: accumulate wk[i] * term into running total with single rounding
		const long double t{ std::fma(wk_[i], term, sum) };

		// Increment actually added to sum
		// Any bits lost in this step are accounted for in the c-update below.
		const long double inc{ t - sum };

		// Neumaier update
		if (std::fabs(sum) >= std::fabs(inc))
			c += (sum - t) + inc;   // lost bits when |sum| ≥ |inc|
		else
			c += (inc - t) + sum;   // lost bits when |inc| > |sum|

		sum = t; // advance running total
	}
	return static_cast<double>(sum + c);
}