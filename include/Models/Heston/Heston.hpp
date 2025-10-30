/**
* Heston.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/Heston/HestonParams.hpp"
#include "Models/Heston/HestonConfig.hpp"
#include "Math/Quadrature/TanHSinH.hpp"
#include "Core/VolSurface.hpp"
#include <complex>
#include <cstdint> 
#include <memory>

class Heston
{
private:

	//--------------------------------------------------------------------------
	// Private member variables
	//--------------------------------------------------------------------------

	VolSurface volSurf_;
	std::shared_ptr<const TanHSinH> quad_;
	HestonParams params_;
	const HestonConfig config_;

	//--------------------------------------------------------------------------
	// Math
	//--------------------------------------------------------------------------

	// Calculate residues arising from the contour shift
	long double getResidues(long double alpha,
		long double F,
		long double K) const noexcept;

	// Get ITM or OTM damping parameter
	long double getAlpha(long double w) const noexcept;

	// Determine contour shift angle
	long double getPhi(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double w) const noexcept;

	// No branch cut characteristic function 
	// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
	static std::complex<long double> charFunction(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		std::complex<long double> u) noexcept;

public:

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	Heston() = delete;

	// Calibration occurs when object is initialized
	explicit Heston(const VolSurface& volSurf, 
		std::shared_ptr<const TanHSinH> quad,
		const HestonConfig& config);

	//--------------------------------------------------------------------------
	// Pricing
	//--------------------------------------------------------------------------

	// European Call Price
	double callPrice(long double kappa,
		long double theta,
		long double sigma,
		long double rho,
		long double v0,
		long double T,
		long double F,
		long double r,
		long double K) const noexcept;
};