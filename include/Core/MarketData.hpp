/**
* MarketData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

struct MarketData
{
	double r; // continuously compounded annualized risk-free rate
	double q; // continuously compounded annualized dividend yield
	double S; // spot price
};