/**
* Heston.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Models/Heston/HestonParams.hpp"
#include "Core/VolSurface.hpp"

class Heston
{
private:

	VolSurface calVolSurf_;
	HestonParams params_;

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	Heston() = delete;

	// Calibration occurs when object is initialized
	Heston(const VolSurface& volSurf);

};