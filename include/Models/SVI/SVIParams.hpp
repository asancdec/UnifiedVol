/**
* SVIParams.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/21/2025
*/

#ifndef SVI_PARAMS_HPP
#define SVI_PARAMS_HPP

// Raw SVI parameters to calibrate per maturity slice
struct SVIParams
{	
	double a;
	double b;
	double rho;
	double m;
	double sigma;
};

#endif // SVI_PARAMS_HPP
