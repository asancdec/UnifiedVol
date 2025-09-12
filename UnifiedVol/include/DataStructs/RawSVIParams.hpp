/**
* RawSVIParams.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/

#ifndef RAW_SVI_PARAMS_HPP
#define RAW_SVI_PARAMS_HPP

// Raw SVI parameters to calibrate per maturity slice
struct RawSVIParams
{
	double a;
	double b;
	double rho;
	double m;
	double sigma;
};

#endif // RAW_SVI_PARAMS_HPP
