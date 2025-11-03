/**
* HestonCalibrator.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/Heston/HestonCalibrator/HestonCalibrator.hpp"

#include <array>

namespace uv
{
	::std::array<double, 5> HestonCalibrator::initGuess() noexcept
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

	::std::array<double, 5> HestonCalibrator::lowerBounds() noexcept
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

	::std::array<double, 5> HestonCalibrator::upperBounds() noexcept
	{
		return
		{
			10.0,   // kappa
			0.5,    // theta
			10.0,   // sigma
			0.999,  // rho
			0.5     // vo
		};
	}
}

