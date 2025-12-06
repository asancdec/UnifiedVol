/**
* Config.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

#include <limits>

namespace uv::models::heston
{
	struct Config
	{
		Real alphaItm;                                    // Damping parameter value when ln(F/K) >= 0
		Real alphaOtm;                                    // Damping parameter value when ln(F/K) < 0 
		Real eps{ std::numeric_limits<Real>::epsilon() };  // Machine Epsilon
	};
}
