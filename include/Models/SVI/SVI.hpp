/**
* SVI.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"
#include "Math/Calibration/NLopt/CalibratorNLopt.hpp"
#include "Models/SVI/SVISlice.hpp"

#include <nlopt.hpp> 
#include <array>      
#include <vector>
#include <tuple>

namespace uv::models::svi
{
	//--------------------------------------------------------------------------
	// Calibration
	//--------------------------------------------------------------------------

	// Main calibration function
	template <nlopt::algorithm Algo>
	std::tuple<std::vector<SVISlice>, VolSurface> calibrate(const VolSurface& mktVolSurf,
		const CalibratorNLopt<5, Algo>& prototype,
		bool isValidateResults = true);

	// g(k) standard function
	double gk(double a, double b, double rho, double m, double sigma, double k) noexcept;

	namespace detail
	{
		//--------------------------------------------------------------------------
		// Forward declarations
		//--------------------------------------------------------------------------

		struct ConstraintCtx;
		struct ObjCtx;
		struct GKPrecomp;

		//--------------------------------------------------------------------------
		// Initial guess and bounds
		//--------------------------------------------------------------------------

		// Initial guess
		std::array<double, 5> initGuess(const SliceData& slice) noexcept;

		// Lower bounds
		std::array<double, 5> lowerBounds(const SliceData& slice) noexcept;

		// Upper bounds
		std::array<double, 5> upperBounds(const SliceData& slice) noexcept;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------

		// Define the minimum total variance constraint: wMin ≥ 0.0 
		template <nlopt::algorithm Algo>
		void addWMinConstraint(CalibratorNLopt<5, Algo>& calibrator) noexcept;

		// Add Roger Lee left wing and right wing max slope constraints
		template <nlopt::algorithm Algo>
		void addMaxSlopeConstraint(CalibratorNLopt<5, Algo>& calibrator) noexcept;

		// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
		template <nlopt::algorithm Algo>
		void addCalendarConstraint(CalibratorNLopt<5, Algo>& calibrator, std::vector<ConstraintCtx>& contexts) noexcept;

		// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
		template <nlopt::algorithm Algo>
		void addConvexityConstraint(CalibratorNLopt<5, Algo>& calibrator, std::vector<double>& kStorage) noexcept;

		// Add objective function with analytical gradient
		template <nlopt::algorithm Algo>
		void setMinObjective(CalibratorNLopt<5, Algo>& calibrator, const ObjCtx& obj) noexcept;

		//--------------------------------------------------------------------------
		// Math functions
		//--------------------------------------------------------------------------

		// SVI total variance for a given log-forward moneyness
		double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

		// g(k) optimized function for calibration
		double gk(const GKPrecomp& p) noexcept;

		// ∇g(k)
		std::array<double, 5> gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

		//--------------------------------------------------------------------------
		// Testing
		//--------------------------------------------------------------------------

		// Check calibration results
		template <nlopt::algorithm Algo>
		void evalCal(const SVISlice& sviSlice,
			const CalibratorNLopt<5, Algo>& calibrator,
			const std::vector<double>& kSlice,
			const std::vector<double>& wKPrevSlice) noexcept;

		// Make a wK slice
		std::vector<double> makewKSlice(const std::vector<double>& kSlice,
			double a, double b, double rho, double m, double sigma) noexcept;
	} // namespace detail
} 

#include "SVI.inl"