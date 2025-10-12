/**
* Heston.hpp
* Author: Alvaro Sanchez de Carlos
*/

struct HestonParams
{
	double kappa;  // Mean reversion speed
	double theta;  // Long term variance
	double sigma;  // Volatility of volatility
	double rho;    // Correlation
	double v0;	   // Initial variance
};