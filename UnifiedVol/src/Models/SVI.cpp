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


SVI::SVI(const VolSurface& mktVolSurf) : mktVolSurf_(mktVolSurf)
{   

    for (const auto& slice : mktVolSurf_.slices_)
    {
        calibrateSlice(slice);
    }
}

void SVI::calibrateSlice(const SliceData& slice)
{

    // Define initial guess 
    std::vector<double> initGuess
    {
        0.5 * slice.minWT(),      // a
        0.1,                      // b
        -0.5,                     // rho
        0.1 ,                     // m
        0.1                       // sigma
    };

    // Lower bound constraints
    std::vector<double> lB{
        1e-4,                     // a
        0.001,                    // b            
        -0.9999,                  // rho
        2.0 * slice.minLogFM(),   // m
        0.01                      // sigma
    };

    // Upper bound constraints
    std::vector<double> uB
    {
        slice.maxWT(),            // a
        1.0,                      // b
        0.9999,                   // rho
        2.0 * slice.maxLogFM(),   // m
        1.0                       // sigma
    };

    // Safety check: clamp initial guesses within lower and upper bounds
    for (size_t i = 0; i < initGuess.size(); ++i) 
    {
        initGuess[i] = std::min(std::max(initGuess[i], lB[i]), uB[i]);
    }

    // NLopt setup (SLSQP)
    nlopt::opt opt(nlopt::LD_SLSQP, 5);
    opt.set_lower_bounds(lB);
    opt.set_upper_bounds(uB);

    // For each log-forward moneyness point, add one convexity constraint g(k) ≥ 0
    for (size_t i = 0; i < slice.logFM_.size(); ++i)
    {   
        // Add one inequality constraint for each point
        opt.add_inequality_constraint
        (
            // Non-capturing lambda. The unary + converts it to a C function pointer (NLopt requirement).
            +[](                  
                unsigned n,       // dimension of x (here 5). Required by NLopt; unused
                const double* x,  // pointer to current parameter vector: [a,b,rho,m,sigma]
                double* grad,     // optional gradient buffer. Present in signature; unused here.
                void* data        // opaque pointer carrying k
             ) -> double 
            {
                (void)n;                  // silence “unused” warning
                (void)grad;               // silence “unused” warning

                // Recover k (log-forward moneyness) from the opaque pointer 
                const double k{ *static_cast<const double*>(data)};

                // NLopt enforces constraint c(x) ≤ 0 -> g(k) ≥ eps
                constexpr double eps{ 1e-6 };
                return eps - SVI::g_k(
                    x[0],   // a
                    x[1],   // b
                    x[2],   // rho
                    x[3],   // m
                    x[4],   // sigma
                    k       // k
                );
            },
            (void*)&slice.logFM_[i],   // opaque address of the current k value (passed to "data" parameter)
            0.0                        // tolerance level
        );
    }

    // Objective data
    struct Obj { const double* Ks; const double* w; size_t n; };
    const auto& Ks = slice.logFM_;
    const auto& wMkt = slice.wT_;
    Obj obj{ Ks.data(), wMkt.data(), Ks.size() };



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

