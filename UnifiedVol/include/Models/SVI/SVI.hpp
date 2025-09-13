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

	// SliceView data
	struct SliceView;

	// No calendar arbitrage constraint context data
	struct ConstraintCtx;

	// g(k) related precomputation data
	struct GKPrecomp;


	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Single slice calibration
	void calibrateSlice(const SliceData& slice, std::vector<double>& wKPrevSlice);

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

	// Add the minimum total variance constraint: wMin ≥ 0.0 
	static void wMinConstraint(nlopt::opt& opt);

	// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
	static void addCalendarNoArbConstraints(nlopt::opt& opt, const SliceView& cal, const double* prevWk);

	// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
	static void addConvexityConstraints(nlopt::opt& opt,const std::vector<double>& kSlices);

	// Add objective function with analytical gradient
	static void objectiveFunction(nlopt::opt& opt, const SliceView& obj);

	//--------------------------------------------------------------------------
	// Math functions
	//--------------------------------------------------------------------------

	// SVI total variance for a given log-forward moneyness
	static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k)
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	// ∇g(k)
	static std::array<double, 5> gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;


public:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	
	
	std::vector<SVISlice> sviSlices_;  // Calibration data
	VolSurface mktVolSurf_             // Market volatility surface

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Delete default constructor
	SVI() = delete;

	// Calibration occurs when object is initialized
	SVI(const VolSurface& mktVolSurf);

	//--------------------------------------------------------------------------
	// Testing
	//--------------------------------------------------------------------------

	void calibrationEvaluation() const;

};


#endif // SVI_HPP