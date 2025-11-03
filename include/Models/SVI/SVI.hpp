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

namespace uv

{
	class SVI
	{
	private:

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
		static ::std::array<double, 5> initGuess(const SliceData& slice) noexcept;

		// Lower bounds
		static ::std::array<double, 5> lowerBounds(const SliceData& slice) noexcept;

		// Upper bounds
		static ::std::array<double, 5> upperBounds(const SliceData& slice) noexcept;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		// Define the minimum total variance constraint: wMin ≥ 0.0 
		template <::nlopt::algorithm Algo>
		static void addWMinConstraint(CalibratorNLopt<5, Algo>& calibrator) noexcept;

		// Add Roger Lee left wing and right wing max slope constraints
		template <::nlopt::algorithm Algo>
		static void addMaxSlopeConstraint(CalibratorNLopt<5, Algo>& calibrator) noexcept;

		// Add no calendar spread arbitrage constraint: Wk_current ≥ Wk_previous
		template <::nlopt::algorithm Algo>
		static void addCalendarConstraint(CalibratorNLopt<5, Algo>& calibrator, ::std::vector<ConstraintCtx>& contexts) noexcept;

		// Add no butterfly arbitrage constraints g(k) ≥ 0.0 
		template <::nlopt::algorithm Algo>
		static void addConvexityConstraint(CalibratorNLopt<5, Algo>& calibrator, ::std::vector<double>& kStorage) noexcept;

		// Add objective function with analytical gradient
		template <::nlopt::algorithm Algo>
		static void setMinObjective(CalibratorNLopt<5, Algo>& calibrator, const ObjCtx& obj) noexcept;

		//--------------------------------------------------------------------------
		// Math functions
		//--------------------------------------------------------------------------
		// SVI total variance for a given log-forward moneyness
		static double wk(double a, double b, double rho, double m, double sigma, double k) noexcept;

		// g(k) optimized function for calibration
		static double gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept; 

		// ∇g(k)
		static ::std::array<double, 5> gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept;

		//--------------------------------------------------------------------------
		// Testing
		//--------------------------------------------------------------------------
		// Check calibration results
		template <::nlopt::algorithm Algo>
		static void evalCal(const SVISlice& sviSlice,
			const CalibratorNLopt<5, Algo>& calibrator,
			const ::std::vector<double>& kSlice,
			const ::std::vector<double>& wKPrevSlice) noexcept;

		// Make a wK slice
		static ::std::vector<double> makewKSlice(const ::std::vector<double>& kSlice,
			double a, double b, double rho, double m, double sigma) noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		SVI() = delete;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		template <::nlopt::algorithm Algo>
		static ::std::tuple<::std::vector<SVISlice>, VolSurface> calibrate(const VolSurface& mktVolSurf,
			const CalibratorNLopt<5, Algo>& prototype,
			bool isValidateResults = true);
	};
}

#include "SVI.inl"