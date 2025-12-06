/**
* MarketData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

namespace uv::core
{
	struct MarketData
	{
		double r; // continuously compounded risk-free rate
		double q; // continuously compounded dividend yield
		double S; // spot price
	};
}