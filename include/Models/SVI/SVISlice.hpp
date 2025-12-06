/**
* SVISlice.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

namespace uv::models::svi
{
	struct SVISlice
	{
		Real T;
		Real a;
		Real b;
		Real rho;
		Real m;
		Real sigma;
	};
}