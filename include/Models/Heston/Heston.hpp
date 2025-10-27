///**
//* Heston.hpp
//* Author: Alvaro Sanchez de Carlos
//*/
//
//#pragma once
//
//#include "Models/Heston/HestonParams.hpp"
//#include "Math/Quadrature/GaussLaguerre.hpp"
//#include "Core/VolSurface.hpp"
//#include <complex>
//#include <cstdint> 
//
//class Heston
//{
//private:
//
//	VolSurface volSurf_;
//	GaussLaguerre glQuadrature_;
//	HestonParams params_;
//
//public:
//
//	//--------------------------------------------------------------------------
//	// Initialization
//	//--------------------------------------------------------------------------
//
//	Heston() = delete;
//
//	// Calibration occurs when object is initialized
//	Heston(const VolSurface& volSurf, const GaussLaguerre& glQuadrature);
//
//	//--------------------------------------------------------------------------
//	// Pricing
//	//--------------------------------------------------------------------------
//
//	// European Call Price
//	double callPrice(double kappa,
//		double theta,
//		double sigma,
//		double rho,
//		double v0,
//		double T,
//		double S,
//		double r,
//		double q,
//		double K) const noexcept;
//
//
//	// No branch cut ( Little Heston Trap) and has analytical gradient (Cui)
//	static std::complex<double> charFunction(double kappa,
//		double theta,
//		double sigma,
//		double rho,
//		double v0,
//		double T,
//		double S,
//		double logSPlusDrift,
//		std::complex<double> u) noexcept;
//
//
//};