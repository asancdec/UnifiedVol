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

struct SVI::SliceView
{
    const double* k;        // Pointer to log-forward moneyness 
    const double* wK;       // Pointer to total variance 
    size_t n;               // Number of data points
    const double* T;        // Pointer to maturity
    const double* vegaMkt;  // Pointer to market Vega

    // Constructor with optional T and vegaMkt
    SliceView(
        const double* k,
        const double* wK,
        size_t n,
        const double* T = nullptr,
        const double* vegaMkt = nullptr)
        : k(k), wK(wK), n(n), T(T), vegaMkt(vegaMkt) {}
};

struct SVI::ConstraintCtx
{
    double k;       // Log-forward moneyness
    double prevWk;  // Total variance of the previous slice
};

struct SVI::GKPrecomp
{
    double x;             // x := k-m
    double R;             // R:= sqrt(x^2 + sigma^2)
    double invR;          // invR := 1 / R
    double wk;            // w(k) := a + b*(rho*x + R)
    double wkD1;          // w'(k) := b * (rho + x/R)
    double wkD1Squared;   // w'(k)^2
    double invRCubed;     // 1/(R^3)
    double sigmaSquared;  // sigma^2
    double wkD2;          // w''(k) := b * sigma^2 / R^3
    double A;             // A := 1 - k * w'/(2 * w)                        
    double B;             // B := 1/w(k) + 1/4

    GKPrecomp(double a, double b, double rho, double m, double sigma, double k) noexcept
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
};


SVI::SVI(const VolSurface& mktVolSurf) : mktVolSurf_(mktVolSurf)
{   
    // Initialize slices vector
    sviSlices_.reserve(mktVolSurf.slices().size());

    // Initialize the first total variance vector
    std::vector<double> wkSlice(mktVolSurf.numStrikes(), 0.0);

    // Calibrate each slice
    for (const auto& slice : mktVolSurf.slices())
    {
        calibrateSlice(slice, wkSlice);

        std::cout << std::fixed << std::setprecision(10) << "wkSlice:";
        for (const double w : wkSlice) std::cout << ' ' << w;
        std::cout << '\n';
    }
}

void SVI::calibrateSlice(const SliceData& slice, std::vector<double>& wKPrevSlice)
{
    // Constant references to slice data
    const std::vector<double>& kSlice{ slice.logFM() };
    const std::vector<double>& wTSlice{ slice.wT() };
    double T{ slice.T() };
    const std::vector<double>& vegaMkt{ slice.vega() };

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

    // Enforce positive minimum total variance constraint
    wMinConstraint(opt);

    // Initialize no calendar arbitrage data instance
    SliceView cal{ kSlice.data(), wKPrevSlice.data(),  wKPrevSlice.size()};

    // Enforce calendar spread arbitrage constraints: Wk_current ≥ Wk_previous
    addCalendarNoArbConstraints(opt, cal, wKPrevSlice.data());

    // Enforce convexity constraints: g(k) ≥ 0
    addConvexityConstraints(opt, kSlice);

    // Initialize objective data instance
    SliceView obj{ kSlice.data(), wTSlice.data(),  wTSlice.size(), &T, vegaMkt.data()};

    // Define objective function with analytical gradient
    objectiveFunction(opt, obj);

    // Stopping criteria
    opt.set_ftol_rel(1e-10);  // Stop when objective stops improving
    opt.set_maxeval(7000);    // Do not evaluate objective/constraints more than 7000 times

    // Solve the problem
    std::vector<double> x{ iG };
    double SSE{ 0.0 };
    nlopt::result res{ opt.optimize(x, SSE) };

    // Extract variables
    const double a{ x[0] };
    const double b{ x[1] };
    const double rho{ x[2] };
    const double m{ x[3] };
    const double sigma{ x[4] };

    // Save calibration results
    sviSlices_.emplace_back(
        SVISlice
        {
        slice.T(),                 
        SVIParams{ a, b, rho, m, sigma }
        });

    // Modify slice total variance vector
    // Purpose: enforce no-arbitrage calendar spread constraint on the next slices 
    std::transform(kSlice.begin(), kSlice.end(), wKPrevSlice.begin(),
        [a, b, rho, m, sigma](double k) noexcept 
        {
            return wk(a, b, rho, m, sigma, k);
        });
}


std::array<double, 5> SVI::initGuess(const SliceData& slice) noexcept
{   
    const double b{0.1};
    const double rho{ -0.5 };
    const double m{ 0.1 };
    const double sigma{ 0.1 };
    const double a{ slice.atmWT() - b * (-rho * m + std::hypot(m, sigma)) };

    return { a, b, rho, m, sigma };
}

std::array<double, 5> SVI::lowerBounds(const SliceData& slice) noexcept
{
    return     
    {
        -2.0,                     // a
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
    3.0                       // sigma
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

void SVI::wMinConstraint(nlopt::opt& opt)
{
    opt.add_inequality_constraint
    ( 
        +[](unsigned n, const double* x, double* grad, void*
            ) -> double
        {
            // Extract variables
            const double a{ x[0] };
            const double b{ x[1] };
            const double rho{ x[2] };
            const double sigma{ x[4] };

            // Precompute variables
            const double S{ std::sqrt(1.0 - rho * rho) };

            // Calculate minimum variance
            const double wMin{ std::fma(b, sigma * S, a) };  // wMin := a + b * sigma * sqrt(1 - rho^2)

            if (grad)
            {
                // ∇c(x) = -∇wMin
                grad[0] = -1.0;                   // -∂wMin/∂a
                grad[1] = -sigma * S;             // -∂wMin/∂b
                grad[2] = b * sigma * (rho / S);  // -∂wMin/∂rho
                grad[3] = 0.0;                    // -∂wMin/∂m
                grad[4] = -b * S;                 // -∂wMin/∂sigma
            }

            // Enforce wMin ≥ eps <-> c(x) = eps - wMin ≤ 0
            return 1e-10 - wMin;
        },
        nullptr,
        1e-10
    );
}

void SVI::addCalendarNoArbConstraints(nlopt::opt& opt, const SliceView& cal, const double* prevWk)
{
    // Store contexts in a static vector so addresses remain valid
    static std::vector<ConstraintCtx> contexts;
    contexts.resize(cal.n);                      
    for (size_t i = 0; i < cal.n; ++i)
    {
        contexts[i] = ConstraintCtx{ cal.k[i], prevWk[i] };
    }

    // For each k, add one no calendar arbitrage constraint
    for (auto& ctx : contexts) 
    {
        opt.add_inequality_constraint(
            +[](unsigned, const double* x, double* grad, void* data) -> double 
            {
                // Reinterpret opaque pointer
                const ConstraintCtx& c{ *static_cast<const ConstraintCtx*>(data) };

                // Extract variables
                const double a{ x[0] };
                const double b{ x[1] };
                const double rho{ x[2] };
                const double m{ x[3] };
                const double sigma{ x[4] };

                // Precompute variables
                const double xi{ c.k - m };
                const double R{ std::hypot(xi, sigma) };
                const double invR{ 1.0 / R };

                // Calculate wK
                const double wK{ a + b * (rho * xi + R) };

                if (grad) 
                {
                    grad[0] = -1.0;                     // -∂w/∂a
                    grad[1] = -(rho * xi + R);          // -∂w/∂b
                    grad[2] = -b * xi;                  // -∂w/∂ρ
                    grad[3] = b * (rho + xi * invR);    // -∂w/∂m
                    grad[4] = -b * (sigma * invR);      // -∂w/∂sigma
                }

                // NLopt expects c(x) ≤ 0.
                // Here:  c(x) = prevWk + eps − w(k)
                // Enforces w(k) ≥ prevWk + eps  (no calendar arbitrage, with small tolerance)
                return c.prevWk + 1e-10 - wK;  
            },
            &ctx,
            1e-10
        );
    }
}

void SVI::addConvexityConstraints(nlopt::opt& opt, const std::vector<double>& kSlice) 
{
    // For each k, add one convexity constraint g(k) ≥ 0
    for (size_t i = 0; i < kSlice.size(); ++i)
    {   
        opt.add_inequality_constraint
        (
            +[](unsigned n, const double* x, double* grad, void* data) -> double
            {
                // Reinterpret opaque pointer
                const double k{ *static_cast<const double*>(data) };

                // Extract variables
                const double a{ x[0] };
                const double b{ x[1] };
                const double rho{ x[2] };
                const double m{ x[3] };
                const double sigma{ x[4] };

                // Build g(K) precompute instance for efficiency
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
                // Here c(x) = -g(k), so -g(k) ≤ tol  
                // Enforces g(k) ≥ -tol  (≈ g(k) ≥ 0 with small slack)
                return -g;
            },
            const_cast<double*>(&kSlice[i]),
            1e-10                  
        );
    }
}

void SVI::objectiveFunction(nlopt::opt& opt, const SliceView& obj)
{   
    opt.set_min_objective
    (
        +[](unsigned n, const double* x, double* grad, void* data) -> double
        {
            // Reinterpret opaque pointer
            const SliceView& obj{ *static_cast<const SliceView*>(data) };

            // Extract variables
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
                // Extract variables
                const double k{ obj.k[i] };
                const double T{ *obj.T };

                // Model data
                const double wKModel{SVI::wk(a,b, rho, m, sigma, k)};
                const double volModel{ std::sqrt(wKModel / T) };

                // Market data
                const double wKMarket{ obj.wK[i] };
                const double volMarket{std::sqrt(wKMarket / T) };

                // Calculate vega weight
                const double weight{ 1.0 / obj.vegaMkt[i]};

                // Calculate weighted SSE between model and market volatility
                const double r{ volModel - volMarket };
                SSE += weight *  r * r;

                if (grad)
                {   
                    // Precompute variables
                    const double xi{ k - m };
                    const double R{ std::hypot(xi, sigma) };
                    const double invR{ 1.0 / R };

                    // Chain rule factor: ∂σ/∂θ = 1/(2*sqrt(T*w)) * ∂w/∂θ
                    const double factor{ 0.5 / std::sqrt(T * wKModel) };

                    // Calculate partial derivatives
                    const double dvolda{ factor };                               // ∂σ/∂a := 1/(2*sqrt(T*w)) * ∂w/∂a   = factor * 1
                    const double dvoldb{ factor * (rho * xi + R) };              // ∂σ/∂b := 1/(2*sqrt(T*w)) * ∂w/∂b   = factor * (rho*xi + R)
                    const double dvolrho{ factor * (b * xi) };                   // ∂σ/∂ρ := 1/(2*sqrt(T*w)) * ∂w/∂ρ   = factor * (b*xi)
                    const double dvoldm{ factor * (-b * (rho + xi * invR)) };    // ∂σ/∂m := 1/(2*sqrt(T*w)) * ∂w/∂m   = factor * ( -b*(rho + xi/R) )
                    const double dvoldsigma{ factor * (b * (sigma * invR)) };    // ∂σ/∂σ := 1/(2*sqrt(T*w)) * ∂w/∂σ   = factor * ( b*(σ/R) )

                    // ∇f += 2 * weight * r * ∂σ/∂θ 
                    g[0] += 2.0 * weight * r * dvolda;
                    g[1] += 2.0 * weight * r * dvoldb;
                    g[2] += 2.0 * weight * r * dvolrho;
                    g[3] += 2.0 * weight * r * dvoldm;
                    g[4] += 2.0 * weight * r * dvoldsigma;
                }
            }

            if (grad)
            {               
                // Provide gradient to NLopt
                const auto mCopy = (std::min)(static_cast<std::size_t>(n), g.size());
                std::copy_n(g.data(), mCopy, grad);
            }
            return SSE;
        },
        const_cast<SliceView*>(&obj)
    );
}

double SVI::wk(double a, double b, double rho, double m, double sigma, double k) noexcept
{   
    const double x {k - m};                                     // x := k-m
    return std::fma(b, (rho * x + std::hypot(x, sigma)), a);    // w(k) := a + b*(rho*x + sqrt(x^2 + sigma^2)) 
}   

double SVI::gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
{   
    return (p.A * p.A) - 0.25 * p.wkD1Squared* p.B + p.wkD2 / 2.0;  // g(k) := A^2 - B * (w')^2 / 4 + w''/2
}

std::array<double, 5> SVI::gkGradient(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
{
    // Precompute variables
    double invR5{ p.invRCubed * p.invR * p.invR };    // 1/R^5
    double wkSquaredInv{ 1.0 / (p.wk * p.wk)};        // 1/w^2
    
    // Calculate partial derivatives of g(k) with respect to w, w', and w''
    double dgdw{ p.A * k * p.wkD1 * wkSquaredInv + 0.25 * p.wkD1Squared * wkSquaredInv }; // ∂g/∂w := A*k*w'/w^2 + (w')^2 / (4 * w^2)
    double dgdw1{ -p.A * k / p.wk - 0.5 * p.wkD1 * p.B};                                  // ∂g/∂w' := -A*k/w - w'*B/2
    double dgdw2{ 0.5 };                                                                  // ∂g/∂w'' := 1/2

    // ---------- ∂w/∂θ ---------
    std::array<double, 5> dw
    {
        1.0,                         // ∂w/∂a   := 1
        rho * p.x + p.R,             // ∂w/∂b   := ρ (k-m) + R
        b * p.x,                     // ∂w/∂ρ   := b (k-m)
        -b * (rho + p.x * p.invR),   // ∂w/∂m   := -b ( ρ + (k-m)/R )
        b * sigma * p.invR           // ∂w/∂σ   := b (σ / R)
    };

    // ---------- ∂w′/∂θ ------
    std::array<double, 5> dw1
    {
        0.0,                                // ∂w′/∂a   := 0
        rho + p.x * p.invR,                 // ∂w′/∂b   := ρ + (k-m)/R
        b,                                  // ∂w′/∂ρ   := b
        -b * p.sigmaSquared * p.invRCubed,  // ∂w′/∂m   := -b * σ^2 / R^3
        -b * p.x * sigma * p.invRCubed      // ∂w′/∂σ   := -b * (k-m) * σ / R^3
    };

    // ---------- ∂w″/∂θ ----------
    std::array<double, 5> dw2
    {
        0.0,                                                                 // ∂w″/∂a   := 0
        p.sigmaSquared * p.invRCubed,                                        // ∂w″/∂b   := σ^2 / R^3
        0.0,                                                                 // ∂w″/∂ρ   := 0
        3.0 * b * p.sigmaSquared * p.x * invR5,                              // ∂w″/∂m   := 3 b σ^2 (k-m) / R^5
        b * (2 * sigma * p.invRCubed - 3.0 * p.sigmaSquared * sigma * invR5) // ∂w″/∂σ := b( 2σ/R^3 - 3σ^3/R^5 )
    };

    // Chain rule: ∇g = (∂g/∂w)∇w + (∂g/∂w1)∇w1 + (∂g/∂w2)∇w2
    std::array<double, 5> dg{};
    for (int j = 0; j < 5; ++j)
    {   
        dg[j] = dgdw * dw[j] + dgdw1 * dw1[j] + dgdw2 * dw2[j];
    }
    return dg;
}