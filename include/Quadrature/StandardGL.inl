/**
* StandardGL.inl
* Author: Alvaro Sanchez de Carlos
*/


template<typename F>
double StandardGL::eval(F&& f) const noexcept
{

	double sum{ 0.0 };   // Running total of the compensated sum
	double c{ 0.0 };     // Compensation accumulator for lost low-order bits

	for (int i = 0; i < N_; ++i)
	{
		// Fused multiply-add: computes sum + wk[i]*f(xk[i]) with a single rounding
		const double t{ std::fma(wk_[i], f(xk_[i]), sum) };

		// Increment actually added to sum
		// Any bits lost in this step are accounted for in the c-update below.
		const double term{ t - sum };

		// Neumaier update
		if (std::fabs(sum) >= std::fabs(term))
			c += (sum - t) + term;   // lost bits when |sum| ≥ |term|
		else
			c += (term - t) + sum;   // lost bits when |term| > |sum|

		sum = t; // advance running total
	}
	return sum + c;
}

