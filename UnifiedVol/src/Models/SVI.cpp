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
#include <array>
#include <iomanip>
#include <limits> 

struct SVI::SliceView
{
    const double* k;        // Pointer to log-forward moneyness 
    const double* wK;       // Pointer to total variance 
    const size_t n;         // Number of data points
};

struct SVI::ConstraintCtx
{
    double k;       // Log-forward moneyness
    double prevWk;  // Total variance of the previous slice
    double eps;     // Epsilon value
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

    explicit GKPrecomp(double a, double b, double rho, double m, double sigma, double k) noexcept
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

SVI::SVI(const VolSurface& mktVolSurf, bool isPrintResults) : config_{}
{   
    // Initialize slices vector
    sviSlices_.reserve(mktVolSurf.slices().size());

    // Initialize the first total variance vector
    std::vector<double> wkSlice(mktVolSurf.numStrikes(), 0.0);

    // Calibrate each slice
    for (const auto& mktSlice : mktVolSurf.slices())
    {
        calibrateSlice(mktSlice, wkSlice, isPrintResults);
    }
}

void SVI::calibrateSlice(const SliceData& mktSlice, std::vector<double>& wKPrevSlice, bool isPrintResults)
{
    // Constant references to slice data
    const std::vector<double>& kSlice{ mktSlice.logFM() };
    const std::vector<double>& wTSlice{ mktSlice.wT() };

	// Size checks
    if (kSlice.size() != wTSlice.size() || kSlice.size() != wKPrevSlice.size())
        throw std::invalid_argument("kSlice/wTSlice/wKPrevSlice size mismatch");

    // Define initial guess 
    std::array<double, 5 > iGArr{ initGuess(mktSlice) };

    // Lower bound constraints
    std::array<double, 5> lBArr{ lowerBounds(mktSlice) };

    // Upper bound constraints
    std::array<double, 5> uBArr{ upperBounds(mktSlice) };

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

    // Enforce Roger Lee wing slope constraints
    addLeeMaxSlopeConstraints(opt);

    // Initialize no calendar arbitrage data instance
    SliceView cal{ kSlice.data(), wKPrevSlice.data(), wKPrevSlice.size()};

    // Vector of constraints context
    std::vector<ConstraintCtx> contexts;
    contexts.resize(cal.n);
    for (size_t i = 0; i < cal.n; ++i)
    {
        contexts[i] = ConstraintCtx{ cal.k[i], wKPrevSlice[i], config_.eps};
    }

    // Enforce calendar spread arbitrage constraints: Wk_current ≥ Wk_previous
    addCalendarNoArbConstraints(opt, contexts);

    // Enforce convexity constraints: g(k) ≥ 0
    addConvexityConstraints(opt, kSlice);

    // Initialize objective data instance
    SliceView obj{ kSlice.data(), wTSlice.data(),  wTSlice.size()};

    // Define objective function with analytical gradient
    objectiveFunction(opt, obj);

    // Stopping criteria
    opt.set_ftol_rel(config_.ftolRel);  // Stop when objective stops improving
    opt.set_maxeval(config_.maxEval);   // Maximum number of evaluations

    // Solve the problem
    std::vector<double> x{ iG };
    double sse{ 0.0 };
    nlopt::result res{ opt.optimize(x, sse) };

    // Extract calibration results
    const double a{ x[0] };
    const double b{ x[1] };
    const double rho{ x[2] };
    const double m{ x[3] };
    const double sigma{ x[4] };
    SVISlice modelSlice{ mktSlice.T(), SVIParams{a, b, rho, m, sigma}};

    // Evaluate calibration
    evalCalib(modelSlice, lBArr, uBArr, sse, mktSlice, kSlice, wKPrevSlice, isPrintResults);

    // Save calibration results
    sviSlices_.emplace_back(modelSlice);

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

std::array<double, 5> SVI::lowerBounds(const SliceData& mktSlice) noexcept
{
    return     
    {
        -10.0,                        // a
        0.001,                       // b            
        -0.9999,                     // rho
        5.0 * mktSlice.minLogFM(),   // m
        0.01                         // sigma
    };
}

std::array<double, 5> SVI::upperBounds(const SliceData& mktSlice) noexcept
{
    return 
    {
    mktSlice.maxWT(),            // a
    2.0,                         // b
    0.9999,                      // rho
    5.0 * mktSlice.maxLogFM(),   // m
    10.0                        // sigma
    };
}

void SVI::clampIG(std::array<double, 5>& iG,
    const std::array<double, 5>& lb,
    const std::array<double, 5> & ub) noexcept
{
    bool clampedAny{ false };

    for (std::size_t i = 0; i < iG.size(); ++i)
    {
        const double before = iG[i];
        const double after = std::clamp(before, lb[i], ub[i]);

        if (after != before)
        {
            clampedAny = true;
            std::cerr << std::fixed << std::setprecision(6)
                << "clampIG: parameter " << paramName(i)
                << " clamped from " << before
                << " to " << after
                << " (lb=" << lb[i] << ", ub=" << ub[i] << ")\n";
            iG[i] = after;
        }
    }
}

void SVI::wMinConstraint(nlopt::opt& opt) const
{   
    const double* epsPtr{ &config_.eps };
    opt.add_inequality_constraint
    ( 
        +[](unsigned, const double* x, double* grad, void*data) -> double
        {
            // Reinterpret opaque pointer
            const double eps {*static_cast<const double*>(data)};

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
            return eps - wMin;
        },
        const_cast<double*>(epsPtr),
        config_.tol
    );
}

void SVI::addLeeMaxSlopeConstraints(nlopt::opt& opt) const
{
    // Right wing: b*(1 + rho) <= 2
    opt.add_inequality_constraint(
        +[](unsigned, const double* x, double* grad, void*) -> double
        {
            // x = [a, b, rho, m, sigma]
            const double b   { x[1] };
            const double rho { x[2] };

            if (grad)
            {
                grad[0] = 0.0;            // d/da
                grad[1] = (1.0 + rho);    // d/db
                grad[2] = b;              // d/drho
                grad[3] = 0.0;            // d/dm
                grad[4] = 0.0;            // d/dsigma
            }

            // c(x) = b*(1 + rho) - 2 <= 0
            return b * (1.0 + rho) - 2.0;
        },
        nullptr,
        config_.tol
    );

    // Left wing: b*(1 - rho) <= 2
    opt.add_inequality_constraint(
        +[](unsigned, const double* x, double* grad, void*) -> double
        {
            const double b   { x[1] };
            const double rho { x[2] };

            if (grad)
            {
                grad[0] = 0.0;            // d/da
                grad[1] = (1.0 - rho);    // d/db
                grad[2] = -b;             // d/drho
                grad[3] = 0.0;            // d/dm
                grad[4] = 0.0;            // d/dsigma
            }

            // c(x) = b*(1 - rho) - 2 <= 0
            return b * (1.0 - rho) - 2.0;
        },
        nullptr,
        config_.tol
    );
}

void SVI::addCalendarNoArbConstraints(nlopt::opt& opt, std::vector<ConstraintCtx>& contexts) const
{
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
                return c.prevWk + c.eps - wK;  
            },
            &ctx,
            config_.tol
        );
    }
}

void SVI::addConvexityConstraints(nlopt::opt& opt, const std::vector<double>& kSlice) const
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
            config_.tol
        );
    }
}

void SVI::objectiveFunction(nlopt::opt& opt, const SliceView& obj) const
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
                const double wKMarket{ obj.wK[i] };

                // Precompute variables
                const double xi{ k - m };
                const double R{ std::hypot(xi, sigma) };
                const double wKModel{ std::fma(b, rho * xi + R, a) };

                // Calculate SSE between model and market variance
                const double r{ wKModel - wKMarket};
                SSE += r * r;

                if (grad)
                {   
                    // Precompute variables
                    const double invR{ 1.0 / R };

                    // Calculate partial derivatives
                    const double dwda{ 1.0 };                    // ∂w/∂a
                    const double dwdb{ (rho * xi + R) };         // ∂w/∂b
                    const double dwdrho{ (b * xi) };             // ∂w/∂ρ
                    const double dwdm{ -b * (rho + xi * invR) }; // ∂w/∂m
                    const double dwdsigma{ b * (sigma * invR) }; // ∂w/∂σ

                    // ∇f += 2 * r * ∂w/∂θ  (variance-SSE)
                    g[0] += 2.0 * r * dwda;
                    g[1] += 2.0 * r * dwdb;
                    g[2] += 2.0 * r * dwdrho;
                    g[3] += 2.0 * r * dwdm;
                    g[4] += 2.0 * r * dwdsigma;
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

void SVI::evalCalib(const SVISlice& modelSlice,
    const std::array<double, 5>& lBArr,
    const std::array<double, 5>& uBArr,
    double sse,
    const SliceData& mktSlice,
    const std::vector<double>& kSlice,
    const std::vector<double>& wKPrevSlice,
    bool isPrintResults) const noexcept
{
    // Extract parameters
    const double T{ modelSlice.T };
    const double a{ modelSlice.sviParams.a };
    const double b{ modelSlice.sviParams.b };
    const double rho{ modelSlice.sviParams.rho };
    const double m{ modelSlice.sviParams.m };
    const double sig{ modelSlice.sviParams.sigma };
    const double SSE{ sse };

    // ---- Print results ----
    if (isPrintResults)
    {
        std::cout << std::fixed << std::setprecision(4)
            << "T=" << T
            << " | a=" << a
            << " b=" << b
            << " rho=" << rho
            << " m=" << m
            << " sigma=" << sig
            << " SSE=" << SSE
            << '\n';
    }

    // ---- Bound touches ----

    // Pack calibrated params for bound checks (uniform init)
    const std::array<double, 5> x{ a, b, rho, m, sig };

    // Near-equality predicate (treat as "hit" if within abs+rel tol)
    auto near = [](double v, double bd) noexcept
        {
            constexpr double abs_eps{ 1e-12 };
            constexpr double rel_eps{ 1e-7 };
            return std::fabs(v - bd) <= (abs_eps + rel_eps * (std::max)(std::fabs(v), std::fabs(bd)));
        };

    for (std::size_t i{ 0 }; i < x.size(); ++i)
    {
        const double v{ x[i] };
        const double lb{ lBArr[i] };
        const double ub{ uBArr[i] };

        if (near(v, lb))
        {
            std::cerr << "WARNING -- " << paramName(i)
                << " hit LOWER bound: v=" << v
                << " (lb=" << lb << ")\n";
        }
        else if (near(v, ub))
        {
            std::cerr << "WARNING -- " << paramName(i)
                << " hit UPPER bound: v=" << v
                << " (ub=" << ub << ")\n";
        }
    }

    // ---- wMin constraint violation ----

    const double S{ std::sqrt((std::max)(0.0, 1.0 - rho * rho)) };
    const double wMin{ std::fma(b, sig * S, a) };
    if (config_.eps - wMin > config_.tol)
    {
        std::cerr << "WARNING -- No-arbitrage constraint violated -- WMIn violated: wMin=" << wMin
            << " < eps=" << config_.eps << "\n";
    }

    // ---- Lee wing-slope constraint violations ----
    // Right wing: b * (1 + rho) <= 2
    const double leeRight{ b * (1.0 + rho) - 2.0 };
    if (leeRight > config_.tol)
    {
        std::cerr << "WARNING -- No-arbitrage constraint violated -- "
            << "Right wing slope: b*(1+rho)=" << (b * (1.0 + rho))
            << " > 2.0  (b=" << b << ", rho=" << rho << ")\n";
    }

    // Left wing: b * (1 - rho) <= 2
    const double leeLeft{ b * (1.0 - rho) - 2.0 };
    if (leeLeft > config_.tol)
    {
        std::cerr << "WARNING -- No-arbitrage constraint violated -- "
            << "Left wing slope: b*(1-rho)=" << (b * (1.0 - rho))
            << " > 2.0  (b=" << b << ", rho=" << rho << ")\n";
    }

    // ---- Convexity g(k) violation ----
    double gMin{ std::numeric_limits<double>::infinity() };
    double kAtMin{ 0.0 };
    std::size_t nViol{ 0 };

    for (const double k : kSlice)
    {
        const GKPrecomp pre{ a, b, rho, m, sig, k };
        const double g{ SVI::gk(a, b, rho, m, sig, k, pre) };
        if (g < gMin) { gMin = g; kAtMin = k; }
        if (g < -config_.tol) { ++nViol; }
    }

    if (gMin < -config_.tol)
    {
        std::cerr << "WARNING -- No-arbitrage constraint violated -- "
            << "Convexity: min g(k)=" << gMin
            << " at k=" << kAtMin
            << "  (violations " << nViol << "/" << kSlice.size()
            << ", tol=" << config_.tol << ")\n";
    }

    // ---- Calendar no-arb: w_curr(k) >= w_prev(k) + eps on current grid ----
    if (wKPrevSlice.size() == kSlice.size())
    {
        std::size_t nCalViol{ 0 };
        double minMargin{ std::numeric_limits<double>::infinity() };
        double kAtWorst{ 0.0 };

        for (std::size_t i{ 0 }; i < kSlice.size(); ++i)
        {
            const double k{ kSlice[i] };
            const double wPrev{ wKPrevSlice[i] + config_.eps };
            const double wCurr{ SVI::wk(a, b, rho, m, sig, k) };

            const double margin{ wCurr - wPrev };            // should be >= 0
            if (margin < minMargin) { minMargin = margin; kAtWorst = k; }
            if (margin < -config_.tol) { ++nCalViol; }
        }

        if (minMargin < -config_.tol)
        {
            std::cerr << "WARNING -- No-arbitrage constraint violated -- "
                << "Calendar: min[w_curr - (w_prev+eps)]=" << minMargin
                << " at k=" << kAtWorst
                << "  (violations " << nCalViol << "/" << kSlice.size()
                << ", tol=" << config_.tol << ", eps=" << config_.eps << ")\n";
        }
    }
}

const char* SVI::paramName(std::size_t i) noexcept
{
    static constexpr const char* names[5] = { "a", "b", "rho", "m", "sigma" };
    return (i < 5) ? names[i] : "?";
}