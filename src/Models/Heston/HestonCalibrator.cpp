/**
* HestonCalibrator.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/Heston/HestonCalibrator.hpp"
#include <ceres/ceres.h>
#include <array>


namespace uv
{
	HestonReport HestonCalibrator::calibrate(const VolSurface mktVolSurf,
		const HestonPricer& pricer)
	{
		//ceres::Problem problem;

		////// Apply in loop
		////for (int i = 0; i < static_cast<int>(lower.size()); ++i)
		////{
		////	problem.SetParameterLowerBound(p, i, lower[i]);
		////	problem.SetParameterUpperBound(p, i, upper[i]);
		////}

	}

	std::array<double, 5> HestonCalibrator::initGuess() noexcept
	{
		return
		{
			2.5,     // kappa
			0.09,    // theta
			0.60,    // sigma
			-0.50,   // rho
			0.09     // vo
		};
	}

	std::array<double, 5> HestonCalibrator::lowerBounds() noexcept
	{
		return
		{
			0.001,   // kappa
			0.001,   // theta
			0.001,   // sigma
			-0.999,  // rho
			0.001    // vo
		};

	}

	std::array<double, 5> HestonCalibrator::upperBounds() noexcept
	{
		return
		{
			10.0,   // kappa
			0.25,   // theta
			2.0,    // sigma
			0.999,  // rho
			0.25    // vo
		};
	}
}
