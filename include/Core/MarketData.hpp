/**
* MarketData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

namespace uv::core
{
	struct MarketData
	{
		Real r; // continuously compounded risk-free rate
		Real q; // continuously compounded dividend yield
		Real S; // spot price
	};
}