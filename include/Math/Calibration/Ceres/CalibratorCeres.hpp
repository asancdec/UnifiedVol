/**
* CalibratorCeres.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Math/Calibration/CalibratorConfig.hpp"
#include <ceres/ceres.h>
#include <array>

namespace uv
{
	template <::std::size_t N>
	class CalibratorCeres
	{
	private:

		//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------
		CalibratorConfig<N> config_;           // Calibration configuration
		::std::array<double, N> lowerBounds_;  // Lower parameter bounds
		::std::array<double, N> upperBounds_;  // Upper parameter bounds

	public:

		::std::array<double, N> x_;            // Parameter block
		::ceres::Problem problem_;             // Ceres problem instance


		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		CalibratorCeres() = delete;
		explicit CalibratorCeres(const CalibratorConfig<N>& config);

		//--------------------------------------------------------------------------
		// Set Initial Guess and Bounds
		//--------------------------------------------------------------------------	
		void setGuessBounds(const ::std::array<double, N>& initGuess,
			const ::std::array<double, N>& lowerBounds,
			const ::std::array<double, N>& upperBounds) noexcept;
	};
}

#include "CalibratorCeres.inl"