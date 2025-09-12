/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/

#ifndef SVI_HPP
#define SVI_HPP

#include "DataStructs/VolSurface.hpp"
#include "DataStructs/RawSVIParams.hpp"
#include "DataStructs/SliceData.hpp"

#include <map>
#include <vector>

class SVI
{
private:

	// Each maturity maps to a slice of SVI parameters
	std::map<double, RawSVIParams> rawParams;

	// Static objective function 
	static double sviObjective(const std::vector<double>& x,std::vector<double>& grad, void* data);

	// Initial guess
	std::vector<double> initGuess(const VolSurface& impVolSurf, const SliceData& sliceData, const double T,
		size_t i, size_t atmInd, size_t prevInd, size_t nextInd) const;

	// Raw SVI to SVI-JW parameters
	static double b(double w, double c, double p);
	static double rho(double p, double w, double b);
	static double a(double vMin, double T, double b, double sigma, double rho);
	static double m(double v, double vMin, double T, double b, double rho, double alpha);
	static double sigma(double m, double alpha);
	static double beta(double rho, double psi, double w, double b);
	static double alpha(double beta);

	// Build model implied volatility surface
	void buildModelSurface(const VolSurface& impVolSurf, const std::vector<std::vector<double>>& k);

	// SVI variance function
	static double sviVariance(double a, double b, double rho, double m, double sigma, double k);

public:

	// Calibrated modelVolSurf
	VolSurface modelVolSurf;

	// Delete default constructor
	SVI() = delete;

	// Custom constructor to calibrate parameters
	SVI(const VolSurface& impVolSurf, unsigned long int maxEval = 1e6);


};


#endif // SVI_HPP