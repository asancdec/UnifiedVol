/**
* HestonParams.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

namespace uv
{
	struct HestonParams
	{
		long double kappa;  // Mean reversion speed
		long double theta;  // Long term variance
		long double sigma;  // Volatility of variance
		long double rho;    // Correlation
		long double v0;     // Initial variance
	};
}