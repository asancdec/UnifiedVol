/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"
#include "Models/SVI/SVISlice.hpp"
#include "Models/Calibration/Config.hpp"
#include "Models/Calibration/Calibrator.hpp"

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
	VolSurface calVolSurf_;            // Calibrated SVI volatility surface

	//--------------------------------------------------------------------------
	// Forward declarations
	//--------------------------------------------------------------------------

	struct ConstraintCtx;
	struct ObjCtx;
	struct GKPrecomp;

	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Single slice calibration
	std::vector<double> calibrateSlice(const SliceData& mktSlice, 
		const std::vector<double>& wKPrevSlice,
		const Config<5>& sviConfig, 
		bool isValidateResults);

	// Initial guess
	static std::array<double, 5> initGuess(const SliceData& slice) noexcept;
	
	// Lower bounds
	static std::array<double, 5> lowerBounds(const SliceData& slice) noexcept;

	// Upper bounds
	static std::array<double, 5> upperBounds(const SliceData& slice) noexcept;

	// Define the minimum total variance constraint: wMin ≥ 0.0 
	static void addWMinConstraint(Calibrator<5>& calibrator) noexcept;

	// Add Roger Lee left wing and right wing max slope constraints
	static void addMaxSlopeConstraint(Calibrator<5>& calibrator) noexcept;

	// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
	static void addCalendarConstraint(Calibrator<5>& calibrator, std::vector<ConstraintCtx>& contexts) noexcept;

	// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
	static void addConvexityConstraint(Calibrator<5>& calibrator, const std::vector<double>& kSlice) noexcept;

	// Add objective function with analytical gradient
	void setMinObjective(Calibrator<5>& calibrator, const ObjCtx& obj) noexcept;

	//--------------------------------------------------------------------------
	// Math functions
	//--------------------------------------------------------------------------

	// SVI total variance for a given log-forward moneyness
	static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	// g(k) optimized function for calibration
	static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept; 

	// ∇g(k)
	static std::array<double, 5> gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

	//--------------------------------------------------------------------------
	// Testing
	//--------------------------------------------------------------------------

	// Check calibration results
	void evalCal(const SVISlice& calibSlice,
		const Config<5>& config,
		const std::vector<double>& kSlice,
		const std::vector<double>& wKPrevSlice) const noexcept;


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
	SVI(const VolSurface& mktVolSurf, bool isValidateResults = true);

	//--------------------------------------------------------------------------
	// Utilities
	//--------------------------------------------------------------------------

	// Get model volatility surface
	const VolSurface& getVolSurf() const noexcept;
};