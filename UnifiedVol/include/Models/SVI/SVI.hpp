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

	std::vector<SVISlice> sviSlices_;  // Calibration data
	SVIConfig config_;                 // Calibration configuration data

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
	void calibrateSlice(const SliceData& mktSlice, std::vector<double>& wKPrevSlice, bool isPrintResults);

	// Initial guess
	static std::array<double, 5> initGuess(const SliceData& mktSlice) noexcept;
	
	// Lower bounds
	static std::array<double, 5> lowerBounds(const SliceData& mktSlice) noexcept;

	// Upper bounds
	static std::array<double, 5> upperBounds(const SliceData& mktSlice) noexcept;

	// Clamp initial guess within the lower and upper bounds
	static void clampIG(std::array<double, 5>& iG,
		const std::array<double, 5>& lb,
		const std::array<double, 5>& ub) noexcept;

	// Add the minimum total variance constraint: wMin ≥ 0.0 
	void wMinConstraint(nlopt::opt& opt) const;

	// Add Roger Lee left wing and right wing max slope constraints
	void addLeeMaxSlopeConstraints(nlopt::opt& opt) const;

	// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
	void addCalendarNoArbConstraints(nlopt::opt& opt, std::vector<ConstraintCtx>& contexts) const;

	// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
	void addConvexityConstraints(nlopt::opt& opt,const std::vector<double>& kSlice) const;

	// Add objective function with analytical gradient
	void objectiveFunction(nlopt::opt& opt, const SliceView& obj) const;

	//--------------------------------------------------------------------------
	// Math functions
	//--------------------------------------------------------------------------

	// SVI total variance for a given log-forward moneyness
	static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k)
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	// ∇g(k)
	static std::array<double, 5> gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	//--------------------------------------------------------------------------
	// Testing
	//--------------------------------------------------------------------------

	// Check calibration results
	void evalCalib(const SVISlice& modelSlice,
		const std::array<double, 5>& lBArr,
		const std::array<double, 5>& uBArr,
		double sse,
		const SliceData& mktSlice,
		const std::vector<double>& kSlice,
		const std::vector<double>& wKPrevSlice,
		bool isPrintResults) const noexcept;

	// Return parameters name based on index
	static const char* paramName(std::size_t i) noexcept;

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// Delete default constructor
	SVI() = delete;

	// Calibration occurs when object is initialized
	SVI(const VolSurface& mktVolSurf, bool isPrintResults = true);

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------
};


#endif // SVI_HPP