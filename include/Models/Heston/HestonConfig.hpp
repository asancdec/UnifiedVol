/**
* HestonConfig.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <limits>

namespace uv
{
	struct HestonConfig
	{
		long double alphaItm;                                    // Damping parameter value when ln(F/K) >= 0
		long double alphaOtm;                                    // Damping parameter value when ln(F/K) < 0 
		double eps{ ::std::numeric_limits<double>::epsilon() };  // Machine Epsilon
	};
}
