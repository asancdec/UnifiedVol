/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#ifndef SVI_HPP
#define SVI_HPP

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVIParams.hpp"

#include <map>
#include <vector>

class SVI
{
private:

	// Calibration

	// Single slice calibration
	void calibrateSlice(const SliceData& slice);

	// SVI total variance 
	static const double wT(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k) function: must be greater or equal to 0 to satisfy butterfly spread no-arbitrage condition
	static const double g_k(double a, double b, double rho, double m, double sigma, double k) noexcept;


public:

	std::map<double, SVIParams> sviParams_;        // Each maturity maps to a slice of raw SVI parameters
	const VolSurface mktVolSurf_;			       // Calibrated volatility surface

	// Initialization

	// Delete default constructor
	SVI() = delete;

	// Calibration occurs when object is initialized
	SVI(const VolSurface& mktVolSurf);




};


#endif // SVI_HPP