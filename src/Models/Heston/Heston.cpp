///**
//* Heston.cpp
//* Author: Alvaro Sanchez de Carlos
//*/
//
//#include "Models/Heston/Heston.hpp"
//#include <cmath>
//#include <numbers>
//#include <iostream>
//#include <random>
//#include <algorithm>
//#include <cstdint>
//
//
//
//Heston::Heston(const VolSurface& volSurf, const TanSinH& glQuadrature) 
//	: volSurf_(volSurf), glQuadrature_(glQuadrature)
//{
//	std::cout << "Heston started\n";
//}
//
//double Heston::callPrice(double kappa,
//	double theta,
//	double sigma,
//	double rho,
//	double v0,
//	double T,
//	double S,
//	double r,
//	double q,
//	double K) const noexcept
//{
//	// Imaginary unit
//	constexpr std::complex<double> i(0.0, 1.0);
//
//	// Common precomputed terms
//	const double logSPlusDrift{ std::log(S) + (r - q) * T };   // ln(S) + (r−q)T
//	const double SXExpMinusQT{ S * std::exp(-q * T) };         // S e^{-qT}
//	const double KXExpMinusRT{ K * std::exp(-r * T) };         // K e^{-rT}
//	const double logK{ std::log(K) };                          // log-strike
//
//	const double invPi{ 1.0 / std::numbers::pi };
//	// Integrand #1:  Re{ e^{-iu logK}/(iu) * φ(u − i, T) } 
//	auto integrand1 = [=](double u) -> double
//		{
//			const std::complex<double> iu{ i * u };
//			return std::real
//			(
//				std::exp(-iu * logK) / iu *
//				Heston::charFunction(kappa, theta, sigma, rho, v0, T, S,
//					logSPlusDrift, std::complex<double>{ u, -1.0 })
//			);
//		};
//	
//	// Integrand #2:  Re{ e^{-iu logK}/(iu) * φ(u, T) }     
//	auto integrand2 = [=](double u) -> double
//		{
//			const std::complex<double> iu{ i * u };
//			return std::real
//			(
//				std::exp(-iu * logK) / iu *
//				Heston::charFunction(kappa, theta, sigma, rho, v0, T, S,
//					logSPlusDrift, std::complex<double>{ u, 0.0})
//			);
//		};
//
//	
//	// Evaluate both integrals ∫₀^∞ simultaneously via Gauss–Laguerre quadrature
//	const auto [integral1, integral2] = glQuadrature_.evalUnweightedBoth(integrand1, integrand2);
//
//	// Final Heston call price  
//	return SXExpMinusQT * (0.5 + invPi * integral1)
//		- KXExpMinusRT * (0.5 + invPi * integral2);
//}
//
//
//
//
//std::complex<double> Heston::charFunction(double kappa,
//	double theta,
//	double sigma,
//	double rho,
//	double v0,
//	double T,
//	double S,
//	double logSPlusDrift,
//	std::complex<double> u) noexcept
//{	
//	//---------------------------------------------------------------------------
//	// Full and fast calibration of	the Heston stochastic volatility model
//	// Yiran Cuia, Sebastian del Baño Rollin, Guido Germano
//	// Eq. (11a): ξ = κ − ρ σ i u
//	// Eq. (11b): d = sqrt( ξ² + σ² (u² + i u) )
//	// Eq. (15b): A₁ = (u² + i u) sinh(d t / 2)
//	// Eq. (15c): A₂ = d cosh(d t / 2) + ξ sinh(d t / 2)
//	// Eq. (15a): A  = A₁ / A₂
//	// Eq. (17b): D  = log d + (κ − d)t/2 − log( (d+ξ)/2 + (d−ξ)/2 e^{−d t} )
//	// Eq. (18): φ(u,t) = exp{ i u (ln S + (r−q)t)
//	//                         − (t κ θ ρ i u)/σ − v₀ A + (2 κ θ / σ²) D }
//	//---------------------------------------------------------------------------
//
//	constexpr std::complex<double> i(0.0, 1.0);
//	const double kappaTimesTheta{ kappa * theta };
//	const double sigmaSquared{ sigma * sigma };
//
//	// Compute auxiliary complex terms: i*u, u^2, (i*u + u^2)
//	const std::complex<double> iu{ i * u };
//	const std::complex<double> uSquared{ u * u };
//	const std::complex<double> iuPlusUSquared{ iu + uSquared };
//
//	// ξ = κ − ρσ i u   (Eq. 11a)
//	const std::complex<double> xi{ kappa - sigma * rho * iu };
//
//	// d = sqrt( ξ² + σ² (u² + i u) )   (Eq. 11b)
//	const std::complex<double> d{ std::sqrt(xi * xi + sigmaSquared * iuPlusUSquared) };
//	const std::complex<double> dTimesTHalf{ d * (T * 0.5) };
//
//	// Evaluate sinh(d t / 2) and cosh(d t / 2) efficiently from a single exp()
//	const std::complex<double> expDT2{ std::exp(dTimesTHalf) };
//	const std::complex<double> invExpDT2{ std::complex<double>(1.0, 0.0) / expDT2 };
//	const std::complex<double> coshDT2{ (expDT2 + invExpDT2) * 0.5 };
//	const std::complex<double> sinhDT2{ (expDT2 - invExpDT2) * 0.5 };
//
//	// Compute A = A1 / A2  (Eq. 15a–15c)
//	const std::complex<double> A1{ iuPlusUSquared * sinhDT2 };
//	const std::complex<double> A2{ d * coshDT2 + xi * sinhDT2 };
//	const std::complex<double> A{ A1 / A2 };
//
//	// exp(-d t) reused in D (Eq. 17b)
//	const std::complex<double> expMinusDT{ invExpDT2 * invExpDT2 };
//
//	// Compute D = log d + (κ − d)t/2 − log((d+ξ)/2 + (d−ξ)/2 e^{−d t})  (Eq. 17b)
//	const std::complex<double> D{
//		std::log(d) + (kappa - d) * (T * 0.5)
//		- std::log((d + xi) * 0.5 + (d - xi) * 0.5 * expMinusDT)
//	};
//
//	// Final characteristic function φ(u,t)  (Eq. 18)
//	return std::exp(
//		iu * logSPlusDrift
//		- (T * kappaTimesTheta * rho * iu) / sigma
//		- v0 * A
//		+ (2.0 * kappaTimesTheta / sigmaSquared) * D);
//}


//
//double Heston::callPrice(double kappa,
//	double theta,
//	double sigma,
//	double rho,
//	double v0,
//	double T,
//	double S,
//	double r,
//	double q,
//	double K) const noexcept
//{
//	// Imaginary unit
//	constexpr std::complex<double> i(0.0, 1.0);
//
//	// Common precomputed terms
//	const double logSPlusDrift{ std::log(S) + (r - q) * T };   // ln(S) + (r−q)T
//	const double SXExpMinusQT{ S * std::exp(-q * T) };         // S e^{-qT}
//	const double KXExpMinusRT{ K * std::exp(-r * T) };         // K e^{-rT}
//	const double logK{ std::log(K) };                          // log-strike
//	const double invPi{ 1.0 / std::numbers::pi };
//
//	// Integrand #1:  Re{ e^{-iu logK}/(iu) * φ(u − i, T) } 
//	auto integrand1 = [=](double u) -> double
//		{
//			const std::complex<double> iu{ i * u };
//			return std::real
//			(
//				std::exp(-iu * logK) / iu *
//				Heston::charFunction(kappa, theta, sigma, rho, v0, T, S,
//					logSPlusDrift, std::complex<double>{ u, -1.0 })
//			);
//		};
//
//	// Integrand #2:  Re{ e^{-iu logK}/(iu) * φ(u, T) }     
//	auto integrand2 = [=](double u) -> double
//		{
//			const std::complex<double> iu{ i * u };
//			return std::real
//			(
//				std::exp(-iu * logK) / iu *
//				Heston::charFunction(kappa, theta, sigma, rho, v0, T, S,
//					logSPlusDrift, std::complex<double>{ u, 0.0})
//			);
//		};
//
//	// Evaluate both integrals ∫₀^∞ simultaneously via Gauss–Laguerre quadrature
//	const auto [integral1, integral2] = glQuadrature_.evalUnweightedBoth(integrand1, integrand2);
//
//	// Final Heston call price  (Eq. 9)
//	return 0.5 * (SXExpMinusQT - KXExpMinusRT)
//		+ (SXExpMinusQT * invPi) * integral1
//		- (KXExpMinusRT * invPi) * integral2;
//}



//std::complex<double> Heston::charFunction(double kappa,
//	double theta,
//	double sigma,
//	double rho,
//	double v0,
//	double T,
//	double S,
//	double logSPlusDrift,
//	std::complex<double> u) noexcept
//{	
//	//---------------------------------------------------------------------------
//	// Full and fast calibration of	the Heston stochastic volatility model
//	// Yiran Cuia, Sebastian del Baño Rollin, Guido Germano
//	// Eq. (11a): ξ = κ − ρ σ i u
//	// Eq. (11b): d = sqrt( ξ² + σ² (u² + i u) )
//	// Eq. (15b): A₁ = (u² + i u) sinh(d t / 2)
//	// Eq. (15c): A₂ = d cosh(d t / 2) + ξ sinh(d t / 2)
//	// Eq. (15a): A  = A₁ / A₂
//	// Eq. (17b): D  = log d + (κ − d)t/2 − log( (d+ξ)/2 + (d−ξ)/2 e^{−d t} )
//	// Eq. (18): φ(u,t) = exp{ i u (ln S + (r−q)t)
//	//                         − (t κ θ ρ i u)/σ − v₀ A + (2 κ θ / σ²) D }
//	//---------------------------------------------------------------------------
//
//	constexpr std::complex<double> i(0.0, 1.0);
//	const double kappaTimesTheta{ kappa * theta };
//	const double sigmaSquared{ sigma * sigma };
//
//	// Compute auxiliary complex terms: i*u, u^2, (i*u + u^2)
//	const std::complex<double> iu{ i * u };
//	const std::complex<double> uSquared{ u * u };
//	const std::complex<double> iuPlusUSquared{ iu + uSquared };
//
//	// ξ = κ − ρσ i u   (Eq. 11a)
//	const std::complex<double> xi{ kappa - sigma * rho * iu };
//
//	// d = sqrt( ξ² + σ² (u² + i u) )   (Eq. 11b)
//	const std::complex<double> d{ std::sqrt(xi * xi + sigmaSquared * iuPlusUSquared) };
//	const std::complex<double> dTimesTHalf{ d * (T * 0.5) };
//
//	// Evaluate sinh(d t / 2) and cosh(d t / 2) efficiently from a single exp()
//	const std::complex<double> expDT2{ std::exp(dTimesTHalf) };
//	const std::complex<double> invExpDT2{ std::complex<double>(1.0, 0.0) / expDT2 };
//	const std::complex<double> coshDT2{ (expDT2 + invExpDT2) * 0.5 };
//	const std::complex<double> sinhDT2{ (expDT2 - invExpDT2) * 0.5 };
//
//	// Compute A = A1 / A2  (Eq. 15a–15c)
//	const std::complex<double> A1{ iuPlusUSquared * sinhDT2 };
//	const std::complex<double> A2{ d * coshDT2 + xi * sinhDT2 };
//	const std::complex<double> A{ A1 / A2 };
//
//	// exp(-d t) reused in D (Eq. 17b)
//	const std::complex<double> expMinusDT{ invExpDT2 * invExpDT2 };
//
//	// Compute D = log d + (κ − d)t/2 − log((d+ξ)/2 + (d−ξ)/2 e^{−d t})  (Eq. 17b)
//	const std::complex<double> D{
//		std::log(d) + (kappa - d) * (T * 0.5)
//		- std::log((d + xi) * 0.5 + (d - xi) * 0.5 * expMinusDT)
//	};
//
//	// Final characteristic function φ(u,t)  (Eq. 18)
//	return std::exp(
//		iu * logSPlusDrift
//		- (T * kappaTimesTheta * rho * iu) / sigma
//		- v0 * A
//		+ (2.0 * kappaTimesTheta / sigmaSquared) * D);
//}