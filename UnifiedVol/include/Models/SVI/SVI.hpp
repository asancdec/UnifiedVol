/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#ifndef SVI_HPP
#define SVI_HPP

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVIParams.hpp"

#include <nlopt.hpp>

#include <map>
#include <vector>
#include <array> 

class SVI
{
private:

	//--------------------------------------------------------------------------
	// Structs
	//--------------------------------------------------------------------------

	// g(k) related precomputed variables
	struct GKPrecomp
	{
		double x;             // x := k-m
		double R;             // R:= sqrt(x^2 + sigma^2)
		double invR;          // invR := 1 / R
		double wk;            // w(k) := a + b*(rho*x + R)
		double wkD1;          // w'(k) := b * (rho + x/R)
		double wkD1Squared;   // w'(k)^2
		double invRCubed;     // 1/(R^3)
		double sigmaSquared;  // sigma^2
		double wkD2;          // w''(k) := b * sigma^2 / R^3
		double A;             // A := 1 - k * w'/(2 * w)                        
		double B;             // B := 1/w(k) + 1/4

		// Constructor
		GKPrecomp(double a, double b, double rho,
			double m, double sigma, double k) noexcept;
	};

	// Objective function data
	struct Obj
	{
		const double* k;   // Pointer to the first element of the forward log-moneyness vector
		const double* wT;  // Pointer to the first element of the total variance vector
		size_t n;          // Number of data points being fit in the objective
	};

	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Single slice calibration
	void calibrateSlice(const SliceData& slice);

	// Initial guess
	static std::array<double, 5> initGuess(const SliceData& slice) noexcept;
	
	// Lower bounds
	static std::array<double, 5> lowerBounds(const SliceData& slice) noexcept;

	// Upper bounds
	static std::array<double, 5> upperBounds(const SliceData& slice) noexcept;

	// Clamp initial guess within the lower and upper bounds
	static void clampIG(std::array<double, 5>& iG,
		const std::array<double, 5>& lb,
		const std::array<double, 5>& ub) noexcept;

	// Add no butterfly spread arbitrage (convexity) constraints g(k) > 0.0 
	static void addConvexityConstraints(nlopt::opt& opt,const std::vector<double>& kSlice,
		double tol);

	// Add objective function with analytical gradient
	static void objectiveFunction(nlopt::opt& opt, const Obj& obj);   

	//--------------------------------------------------------------------------
	// Math functions
	//--------------------------------------------------------------------------

	// SVI total variance for a given log-forward moneyness
	static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k) function: must be greater or equal to 0 to satisfy butterfly spread no-arbitrage condition
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	// Analytical solution to the g(k) constraint gradient using the chain rule
	static std::array<double, 5>  gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;


public:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	
	
	std::map<double, SVIParams> sviParams_;        // Each maturity maps to a slice of raw SVI parameters
	const VolSurface mktVolSurf_;			       // Calibrated volatility surface

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Delete default constructor
	SVI() = delete;

	// Calibration occurs when object is initialized
	SVI(const VolSurface& mktVolSurf);
};


#endif // SVI_HPP