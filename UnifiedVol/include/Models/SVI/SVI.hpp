/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#ifndef SVI_HPP
#define SVI_HPP

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVISlice.hpp"
#include "Models/SVI/SVIConfig.hpp"

#include <nlopt.hpp>

#include <stdexcept>
#include <map>
#include <vector>
#include <array> 

class SVI
{
private:

	//--------------------------------------------------------------------------
	// Member variables
	//--------------------------------------------------------------------------	

	std::vector<SVISlice> sviSlices_;  // Calibration data and parameters
	SVIConfig config_;                 // Calibration configuration data
	VolSurface calVolSurf_;            // Calibrated SVI volatility surface

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
	std::vector<double> calSlice(const SliceData& mktSlice, const std::vector<double>& wKPrevSlice, bool isPrintResults);

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
	void addWMinConstr(nlopt::opt& opt) const;

	// Add Roger Lee left wing and right wing max slope constraints
	void addMaxSlopeConstr(nlopt::opt& opt) const;

	// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
	void addCalendarConstr(nlopt::opt& opt, std::vector<ConstraintCtx>& contexts) const;

	// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
	void addConvexityConstr(nlopt::opt& opt,const std::vector<double>& kSlice) const;

	// Add objective function with analytical gradient
	void addObjFunc(nlopt::opt& opt, const SliceView& obj) const;

	//--------------------------------------------------------------------------
	// Math functions
	//--------------------------------------------------------------------------

	// SVI total variance for a given log-forward moneyness
	static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k) optimized version for calibration
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept; 

	// ∇g(k)
	static std::array<double, 5> gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	//--------------------------------------------------------------------------
	// Testing
	//--------------------------------------------------------------------------

	// Check calibration results
	void evalCal(const SVISlice& modelSlice,
		const std::array<double, 5>& lBArr,
		const std::array<double, 5>& uBArr,
		double sse,
		const std::vector<double>& kSlice,
		const std::vector<double>& wKPrevSlice,
		bool isPrintResults) const noexcept;

	// Return parameters name based on index
	static const char* paramName(std::size_t i) noexcept;

	// Make a wK slice
	static std::vector<double> makewKSlice(const std::vector<double>& kSlice,
		double a, double b, double rho, double m, double sigma) noexcept;

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Delete default constructor
	SVI() = delete;

	// Calibration occurs when object is initialized
	SVI(const VolSurface& mktVolSurf, bool isPrintResults = true);

	//--------------------------------------------------------------------------
	// Math functions 
	//--------------------------------------------------------------------------
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------

	// Get model volatility surface
	const VolSurface& getVolSurf() const noexcept;
};


#endif // SVI_HPP