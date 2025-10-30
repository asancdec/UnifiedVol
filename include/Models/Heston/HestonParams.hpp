/**
* HestonParams.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

struct HestonParams
{
	long double kappa;  // Mean reversion speed
	long double theta;  // Long term variance
	long double sigma;  // Volatility of volatility
	long double rho;    // Correlation
};