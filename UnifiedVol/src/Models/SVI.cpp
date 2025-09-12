/**
* SVI.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#include "Models/SVI/SVI.hpp"

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>


SVI::GKPrecomp::GKPrecomp(double a, double b, double rho, double m, double sigma, double k) noexcept
{
    this->x = k - m;                                          // x := k-m
    this->R = std::hypot(x, sigma);                           // R:= sqrt(x^2 + sigma^2)
    this->invR = 1.0 / R;                                     // invR := 1 / R
    this->wk = std::fma(b, (rho * x + R), a);                 // w(k) := a + b*(rho*x + R)
    this->wkD1 = b * (rho + x * invR);                        // w'(k) := b * (rho + x/R)
    this->wkD1Squared = wkD1 * wkD1;                          // w'(k)^2
    this->invRCubed = invR * invR * invR;                     // 1/(R^3)
    this->sigmaSquared = sigma * sigma;                       // sigma^2
    this->wkD2 = b * sigmaSquared * invRCubed;                // w''(k) := b * sigma^2 / R^3
    this->A = 1.0 - 0.5 * k * wkD1 / wk;                      // A := 1 - k * w'/(2 * w)                        
    this->B = (1.0 / wk) + 0.25;                              // B := 1/w(k) + 1/4
}


SVI::SVI(const VolSurface& mktVolSurf) : mktVolSurf_(mktVolSurf)
{   

    //for (const auto& slice : mktVolSurf_.slices_)
    //{
    //    calibrateSlice(slice);
    //}

   
    for (const auto& slice : mktVolSurf_.slices()) 
    {
        calibrateSlice(slice);
        break; // only the first one
    }
}

void SVI::calibrateSlice(const SliceData& slice)
{
    // Store slice data as constant references
    const std::vector<double>& kSlice{ slice.logFM() };
    const std::vector<double>& wTSlice{ slice.wT() };

    // Define initial guess 
    std::array<double, 5 > iGArr{ initGuess(slice) };

    // Lower bound constraints
    std::array<double, 5> lBArr{ lowerBounds(slice) };

    // Upper bound constraints
    std::array<double, 5> uBArr{ upperBounds(slice) };

    // Clamp initial guess within lower and upper bounds
    clampIG(iGArr, lBArr, uBArr);

    // Convert arrays into vectors for NlOpt
    std::vector<double> iG(iGArr.begin(), iGArr.end());
    std::vector<double> lB(lBArr.begin(), lBArr.end());
    std::vector<double> uB(uBArr.begin(), uBArr.end());

    // NLopt setup (SLSQP)
    nlopt::opt opt(nlopt::LD_SLSQP, 5);
    opt.set_lower_bounds(lB);
    opt.set_upper_bounds(uB);

    // For each log-forward moneyness point, add one convexity constraint g(k) ≥ 0
    addConvexityConstraints(opt, kSlice, 1e-6);

    // Initialize objective data instance
    Obj obj{ kSlice.data(), wTSlice.data(),  wTSlice.size() };

    // Define objective function with analytical gradient
    objectiveFunction(opt, obj);

    // Stopping criteria
    opt.set_ftol_rel(1e-8); // Stop  when objective stops improving
    opt.set_maxeval(5000);  // Do not evaluate objective/constraints more than 4000 times


    std::vector<double> x{ iG };
    double SSE{0.0};
    try {
        nlopt::result res = opt.optimize(x, SSE);
    }
    catch (const std::exception& e) 
    {
        std::cerr << "NLopt error: " << e.what() << '\n';

        return; // or handle
    }

    std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << '\n';
}


std::array<double, 5> SVI::initGuess(const SliceData& slice) noexcept
{
    return     
    {
        0.5 * slice.minWT(),      // a
        0.1,                      // b
        -0.5,                     // rho
        0.1 ,                     // m
        0.1                       // sigma
    };
}

std::array<double, 5> SVI::lowerBounds(const SliceData& slice) noexcept
{
    return     
    {
        1e-4,                     // a
        0.001,                    // b            
        -0.9999,                  // rho
        2.0 * slice.minLogFM(),   // m
        0.01                      // sigma
    };
}

std::array<double, 5> SVI::upperBounds(const SliceData& slice) noexcept
{
    return 
    {
    slice.maxWT(),            // a
    1.0,                      // b
    0.9999,                   // rho
    2.0 * slice.maxLogFM(),   // m
    1.0                       // sigma
    };
}

void SVI::clampIG(std::array<double, 5>& iG,
    const std::array<double, 5>& lb,
    const std::array<double, 5> & ub) noexcept
{
    for (std::size_t i = 0; i < iG.size(); ++i) 
    {
        iG[i] = std::clamp(iG[i], lb[i], ub[i]);
    }
}

void SVI::addConvexityConstraints(nlopt::opt& opt,
    const std::vector<double>& kSlice,
    double tol) 
{
    // For each log-forward moneyness point, add one convexity constraint g(k) ≥ 0
    for (size_t i = 0; i < kSlice.size(); ++i)
    {   
        // Add one inequality constraint for each point
        opt.add_inequality_constraint
        (
            // Non-capturing lambda. 
            // The unary + converts it to a C function pointer (NLopt requirement).
            +[](
                unsigned n,       // dimension of x (here 5)
                const double* x,  // pointer to current parameter vector: [a,b,rho,m,sigma]
                double* grad,     // gradient buffer
                void* data        // opaque pointer carrying k
                ) -> double
            {
                // Recover parameters
                const double k{ *static_cast<const double*>(data) };
                const double a{ x[0] };
                const double b{ x[1] };
                const double rho{ x[2] };
                const double m{ x[3] };
                const double sigma{ x[4] };

                // Build precompute from the CURRENT x and this k
                const SVI::GKPrecomp p{ a, b, rho, m, sigma, k };
                const double g {SVI::gk(a, b, rho, m, sigma, k, p)};

                if (grad) 
                {
                    // ∇g = [dg/da, dg/db, dg/dρ, dg/dm, dg/dσ]
                    const auto dg{ SVI::gkGradient(a, b, rho, m, sigma, k, p) };
                    for (unsigned j = 0; j < 5 && j < n; ++j)
                    {
                        grad[j] = -dg[j]; // c(x) = -g(k) 
                    }
                }

                // NLopt enforces c(x) ≤ tol. 
                // With c(x) = -g(k) this means: -g(k) ≤ tol <-> g(k) ≥ -tol.
                return -g;
            },
            reinterpret_cast<void*>(const_cast<double*>(&kSlice[i])),
            tol                  
        );
    }
}

void SVI::objectiveFunction(nlopt::opt& opt, const Obj& obj)
{   
    // Set minimization objective (SSE)
    opt.set_min_objective
    (
        // Non-capturing lambda. 
        // The unary + converts it to a C function pointer (NLopt requirement).
        +[](
            unsigned n,       // dimension of x (here 5).
            const double* x,  // pointer to current parameter vector: [a,b,rho,m,sigma]
            double* grad,     // gradient buffer
            void* data        // opaque pointer to Obj { Ks, w, n } 
            ) -> double
        {
            // Recover Obj from the opaque pointer 
            const Obj& obj = *static_cast<const Obj*>(data);

            // Recover parameters
            const double a{ x[0] };
            const double b{ x[1] };
            const double rho{ x[2] };
            const double m{ x[3] };
            const double sigma{ x[4] };

            // Initialize variable to store total SSE
            double SSE{ 0.0 };

            // Accumulate ∇f here
            std::array<double, 5> g{}; 

            for (size_t i = 0; i < obj.n; ++i)
            {
                // Evaluate total variance
                const double wT{SVI::wk(a,b, rho, m, sigma, obj.k[i])};

                // Calculate SSE between model and market wT
                const double r{ wT - obj.wT[i] };
                SSE += r * r;

                // Precompute variables
                const double xi{ obj.k[i] - m };
                const double R{ std::hypot(xi, sigma) };
                const double invR{ 1.0 / R };

                // Accumulate analytical gradient of the objective function 
                if (grad) 
                {
                    std::array<double, 5> dw
                    {
                        1.0,                   // ∂w/∂a
                        rho * xi + R,          // ∂w/∂b
                        b * xi,                // ∂w/∂ρ
                        -b *(rho + xi * invR), // ∂w/∂m
                        b* (sigma * invR)      // ∂w/∂sigma
                    };
                    
                    for (int j = 0; j < 5 && j < static_cast<int>(n); ++j)
                    {   
                        // ∇f += 2 r ∇w
                        g[j] += 2.0 * r * dw[j];
                    }
                }
            }

            // Provide gradient to NLopt
            if (grad)
            {
                const auto mCopy = (std::min)(static_cast<std::size_t>(n), g.size());
                std::copy_n(g.data(), mCopy, grad);
            }

            return SSE;
        },
        const_cast<Obj*>(&obj)   
    );
}

double SVI::wk(double a, double b, double rho, double m, double sigma, double k) noexcept
{   
    // Precompute parameters
    const double x {k - m};                     // x := k-m
    const double R{ std::hypot(x, sigma) };     // R:= sqrt(x^2 + sigma^2)

    // Calculate and return total variance
    return std::fma(b, (rho * x + R), a);       // w(k) := a + b*(rho*x + R) 
}

double SVI::gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
{   
    return (p.A * p.A) - 0.25 * p.wkD1Squared* p.B + p.wkD2 / 2.0;   // g(k) := A^2 - B * (w')^2 / 4 + w''/2
}

std::array<double, 5> SVI::gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
{
    // Precomputations
    double invR5{ p.invRCubed * p.invR * p.invR };    // 1/R^5
    double wkSquaredInv{ 1.0 / (p.wk * p.wk)};       // 1/w^2
    
    // Calculate partial derivatives of g(k) with respect to w, w1, and w2
    double dgdw{ p.A * k * p.wkD1 * wkSquaredInv + 0.25 * p.wkD1Squared * wkSquaredInv }; // ∂g/∂w := A*k*w'/w^2 + (w')^2 / (4 * w^2)
    double dgdw1{ -p.A * k / p.wk - 0.5 * p.wkD1 * p.B};                                  // ∂g/∂w1 := -A*k/w - w'*B/2
    double dgdw2{ 0.5 };                                                                  // ∂g/∂w2 := 1/2

    // Calculate gradient of w
    std::array<double, 5> dw
    {
        1.0,                       // a
        rho * p.x + p.R,           // b
        b * p.x,                   // rho
        -b * (rho + p.x * p.invR), // m
        b * sigma * p.invR         // sigma
    };

    // Calculate gradient of w1
    std::array<double, 5> dw1
    {
        0.0,                                // a
        rho + p.x * p.invR,                 // b
        b,                                  // rho
        -b * p.sigmaSquared * p.invRCubed,  // m
        -b * p.x * sigma * p.invRCubed      // sigma
    };

    // Calculate gradient of w2
    std::array<double, 5> dw2
    {
        0.0,                                                                // a
        p.sigmaSquared * p.invRCubed,                                       // b
        0.0,                                                                // rho
        3.0 * b * p.sigmaSquared * p.x * invR5,                             // m
        b*(2 * sigma * p.invRCubed - 3.0*p.sigmaSquared * sigma * invR5)    // sigma
    };


    // Calculate the gradient of g(k) using the chain rule
    std::array<double, 5> dg{};
    for (int j = 0; j < 5; ++j)
    {   
        // ∇g = (∂g/∂w)∇w + (∂g/∂w1)∇w1 + (∂g/∂w2)∇w2
        dg[j] = dgdw * dw[j] + dgdw1 * dw1[j] + dgdw2 * dw2[j];

    }
    return dg;
}