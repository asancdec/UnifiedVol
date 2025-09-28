/**
* StandardGL.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/27/2025
*/

#ifndef STANDARDGL_HPP
#define STANDARDGL_HPP

#include <vector>
#include <functional> 

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
 
	// Evaluate ∫ f(x) e^{-x} dx ≈ Σ w_i f(x_i)  (α = 0 case)
	double eval(const std::function<double(double)>& f) const;

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------
	void printGrid() const noexcept;

};

#endif // STANDARDGL_HPP