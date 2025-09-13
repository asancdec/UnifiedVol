/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#ifndef SVI_HPP
#define SVI_HPP

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVISlice.hpp"

#include <nlopt.hpp>

#include <map>
#include <vector>
#include <array> 

class SVI
{
private:

	//--------------------------------------------------------------------------
	// Forward declarations
	//--------------------------------------------------------------------------
	struct Obj;
	struct GKPrecomp;

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

	// Add the minimum total variance constraint: WMin ≥ 0.0 
	static void wTMinConstraint(nlopt::opt& opt);

	// Add no butterfly spread arbitrage (convexity) constraints g(k) ≥ 0.0 
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
	static std::array<double, 5> gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;


public:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	
	
	std::vector<SVISlice> sviSlices_;              // Calibrated slice
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