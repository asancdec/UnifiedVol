/**
* Config.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"
#include <array>    
#include <string_view>
#include <cstddef>

namespace uv::cal::nlopt
{
	template <std::size_t N>
	struct Config
	{
		double eps;                                       // Small epsilon used for inequality tolerances
		double tol;                                       // Constraint tolerance passed to NLopt (|c(x)| ≤ tol considered satisfied)
		double ftolRel;                                   // Relative stopping criterion for objective improvement (NLopt ftol_rel)
		unsigned int maxEval;                           // Maximum number of function evaluations allowed during optimization
		std::array<std::string_view, N> paramNames;     // Parameter names
	};
}
