/**
* CalibratorCeres.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Math/Calibration/Ceres/CeresConfig.hpp"
#include "CeresPolicy.hpp"

#include <ceres/ceres.h>
#include <array>

namespace uv
{
	template <std::size_t N, typename Policy = CeresPolicy<>>
	class CalibratorCeres
	{
	private:

		//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------
		CeresConfig<N> config_;                // Calibration configuration
		std::array<double, N> lowerBounds_;  // Lower parameter bounds
		std::array<double, N> upperBounds_;  // Upper parameter bounds
		std::array<double, N> x_;            // Parameter block
		ceres::Problem problem_;             // Ceres problem instance

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		CalibratorCeres() = delete;
		explicit CalibratorCeres(const CeresConfig<N>& config);

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------	
		// Set initial guess and bounds
		void setGuessBounds(const std::array<double, N>& initGuess,
			const std::array<double, N>& lowerBounds,
			const std::array<double, N>& upperBounds) noexcept;
		
		// Set differentiation cost function
		void addAnalyticResidual(std::unique_ptr<ceres::CostFunction> cf) noexcept;

		// Solve the problem
		std::array<double, N> optimize();
	};
}

#include "CalibratorCeres.inl"