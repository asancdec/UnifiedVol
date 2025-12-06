/**
* Params.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

namespace uv::models::heston
{
	struct Params
	{
		Real kappa;  // Mean reversion speed
		Real theta;  // Long term variance
		Real sigma;  // Volatility of variance
		Real rho;    // Correlation
		Real v0;     // Initial variance
	};
}