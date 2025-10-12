/**
* GaussLaguerre.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once 

#include <vector>
#include <cmath> 

class GaussLaguerre
{
private:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	

	const int N_;                  // Number of points in the quadrature
	const long double alpha_;	   // Shape parameter
	std::vector<long double> xk_;  // Nodes vector
	std::vector<long double> wk_;  // Weights vector
	//--------------------------------------------------------------------------
	// Math
	//--------------------------------------------------------------------------

	// Evaluate derivative L'_N(α)(x)
	long double Lprime(const long double xi) const noexcept;

	// Gauss–Laguerre weight at node x_i 
	long double weight(const long double xi) const noexcept;

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Constructor uses Golub and Welsh algorithm to generate nodes
	// Calculates weights using derivative of Laguerre polynomial
	GaussLaguerre(const int N = 64, const double = 0.0);

	//--------------------------------------------------------------------------
	// Math
	//--------------------------------------------------------------------------

	// Evaluation a function's numeric integral using GL-Quadrature
	template<typename F>
	double eval(F&& f) const;

	// Evaluation a function's unweighted numeric integral using GL-Quadrature
	template<typename F>
	double evalUnweighted(F&& f) const;

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------
	void printGrid() const noexcept;

};

#include "GaussLaguerre.inl"

