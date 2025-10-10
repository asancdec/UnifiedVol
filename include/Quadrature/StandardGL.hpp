/**
* StandardGL.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once 

#include <vector>
#include <cmath> 

class StandardGL
{
private:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	

	const int N_;             // Number of points in the quadrature
	std::vector<double> xk_;  // Nodes vector
	std::vector<double> wk_;  // Weights vector

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Constructor uses Golub and Welsh algorithm to generate weights and nodes
	StandardGL(const int N = 64, const double alpha = 0.0);

	//--------------------------------------------------------------------------
	// Functionality
	//--------------------------------------------------------------------------
 
	// Evaluation a function's numeric integral using GL-Quadrature
	// Neumaier compensated summation
	template<typename F>
	double eval(F&& f) const noexcept;

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------
	void printGrid() const noexcept;

};

#include "StandardGL.inl"

