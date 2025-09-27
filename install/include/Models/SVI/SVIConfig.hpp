/**
* SVIConfig.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/21/2025
*/

#ifndef SVI_CONFIG_HPP
#define SVI_CONFIG_HPP

// SVI calibration configuration data
struct SVIConfig
{
	double eps{1e-12};
	double tol{1e-9};
	double ftolRel{ 1e-10 };
	unsigned int maxEval{ 10000 };
};

#endif // SVI_CONFIG_HPP

