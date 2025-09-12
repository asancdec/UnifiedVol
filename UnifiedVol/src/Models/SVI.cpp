/**
* SVI.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#include "Models/SVI/SVI.hpp"

#include <nlopt.hpp>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


SVI::SVI(const VolSurface& mktVolSurf) : modelVolSurf(mktVolSurf)
{
    // Define initial guesses
}

const double SVI::wT(double a, double b, double rho, double m, double sigma, double k) noexcept
{   
    // Precompute parameters
    const double x {k - m};
    const double R{ std::hypot(x, sigma) };

    // Calculate and return total variance
    return std::fma(b, (rho * x + R), a);
}

const double SVI::g_k(double a, double b, double rho, double m, double sigma, double k) noexcept
{   
    // Precompute parameters
    const double x { k - m };                                      // x := k-m
    const double R{ std::hypot(x, sigma) };                        // sqrt(x^2 + sigma^2)
    const double invR{ 1.0 / R };

    // Calculate total variance and its first and second derivative with respect to k
    const double w{ std::fma(b, (rho * x + R), a) };               // w(k) := a + b*(rho*x + R)
    const double w1{ b * (rho + x * invR) };                       // w'(k)
    const double w2{ b * (sigma * sigma) * (invR * invR * invR) }; // w''(k)

    // Precompute g(k) terms
    const double A{ 1.0 - 0.5 * k * (w1 / w) };                    // 1 - k w'/(2w)
    const double B{ (1.0 / w) + 0.25 };                            // 1/w + 1/4

    // Calculate and return g(k) 
    return (A * A) - 0.25 * (w1 * w1) * B + w2;
}

