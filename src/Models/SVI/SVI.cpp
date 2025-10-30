/**
* SVI.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/SVI/SVI.hpp"
#include "Utils/Log.hpp"
#include <nlopt.hpp>
#include <algorithm>
#include <cmath>
#include <limits> 
#include <memory>
#include <format>

struct SVI::ConstraintCtx
{
    double k;       // Log-forward moneyness
    double prevWk;  // Total variance of the previous slice
    double eps;     // Epsilon value
};

struct SVI::ObjCtx 
{
    const double* k;   // pointer to log-forward moneyness
    const double* wK;  // pointer to total variance
    std::size_t   n;   // number of points
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

SVIReport SVI::calibrate(const VolSurface& mktVolSurf,
    const Calibrator<5>& prototype,
    bool isValidateResults)
{   
    // Copy market volatility surface
    VolSurface sviVolSurf{ mktVolSurf };

    // Initialize vectors
    std::vector<SVISlice> sviSlices;
    sviSlices.reserve(sviVolSurf.numSlices());
    std::vector<double> wkSlice(sviVolSurf.numStrikes(), 0.0);

    // Calibrate each slice
    for (auto& slice : sviVolSurf.slices())
    {   
        // Constant references to slice data
        const std::vector<double>& kSlice{ slice.logFM() };
        const std::vector<double>& wKSlice{ slice.wT() };

        // Initialize calibrator instance
        Calibrator<5> calibrator{ prototype.fresh() };

        // Set initial guess and bounds
        calibrator.setGuessBounds(
            initGuess(slice),
            lowerBounds(slice),
            upperBounds(slice)
        );

        // Enforce positive minimum total variance constraint
        addWMinConstraint(calibrator);

        // Enforce Roger Lee wing slope constraints
        addMaxSlopeConstraint(calibrator);

        // Vector of constraints context
        std::vector<ConstraintCtx> contexts;
        contexts.resize(kSlice.size());
        for (std::size_t i = 0; i < kSlice.size(); ++i)
        {
            contexts[i] = ConstraintCtx
            {
                kSlice[i],
                wkSlice[i],
                calibrator.eps()
            };
        }

        // Enforce calendar spread arbitrage constraints: Wk_current ≥ Wk_previous
        addCalendarConstraint(calibrator, contexts);

        // Enforce convexity constraints: g(k) ≥ 0
        std::vector<double> kStorage{ kSlice };
        addConvexityConstraint(calibrator, kStorage);

        // Objective function contexts
        ObjCtx obj{ kSlice.data(), wKSlice.data(), kSlice.size() };

        // Define objective function with analytical gradient
        setMinObjective(calibrator, obj);

        // Solve the optimization problem
        std::vector<double> params{ calibrator.optimize() };

        // Extract calibration results
        double T{ slice.T() };
        double a{ params[0] };
        double b{ params[1] };
        double rho{ params[2] };
        double m{ params[3] };
        double sigma{ params[4] };

        SVISlice sviSlice{ T, a, b, rho, m, sigma };

        // Evaluate calibration
        if (isValidateResults) evalCal(sviSlice, calibrator, kSlice, wkSlice);

        // Save calibration parameter results
        sviSlices.emplace_back(std::move(sviSlice));

        // Update calculated variances use them on the next slice calibration
        wkSlice = makewKSlice(kSlice, a, b, rho, m, sigma);

        // Set the slice of the calibrated volatility surface
        slice.setWT(wkSlice);
    }
    return { std::move(sviSlices), std::move(sviVolSurf) };
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
        -10.0,                       // a
        0.001,                       // b            
        -0.9999,                     // rho
        5.0 * slice.minLogFM(),      // m
        0.01                         // sigma
    };
}

std::array<double, 5> SVI::upperBounds(const SliceData& slice) noexcept
{
    return 
    {
    slice.maxWT(),              // a
    2.0,                        // b
    0.9999,                     // rho
    5.0 * slice.maxLogFM(),     // m
    10.0                        // sigma
    };
}

void SVI::addWMinConstraint(Calibrator<5>& calibrator) noexcept
{
    const double* epsPtr{ &calibrator.eps() };

    calibrator.addInequalityConstraint(
        +[](unsigned /*n*/, const double* x, double* grad, void* data) -> double
        {
            // Reinterpret opaque pointer
            const double eps{ *static_cast<const double*>(data) };

            // Extract variables
            const double a = x[0];
            const double b = x[1];
            const double rho = x[2];
            const double sigma = x[4];

            // Precompute
            const double S = std::sqrt(1.0 - rho * rho);

            // w_min = a + b * sigma * sqrt(1 - rho^2)
            const double wMin = std::fma(b, sigma * S, a);

            if (grad) 
            {
                // c(x) = eps - wMin  =>  ∇c = -∇wMin
                grad[0] = -1.0;                   // -∂wMin/∂a
                grad[1] = -sigma * S;             // -∂wMin/∂b
                grad[2] = b * sigma * (rho / S);  // -∂wMin/∂rho (since ∂S/∂rho = -rho/S)
                grad[3] = 0.0;                    // -∂wMin/∂m
                grad[4] = -b * S;                 // -∂wMin/∂sigma
            }

            // Enforce wMin ≥ eps  ⇔  c(x) = eps - wMin ≤ 0
            return eps - wMin;
        },
        const_cast<double*>(epsPtr) // nlopt API takes void*
    );
}

void SVI::addMaxSlopeConstraint(Calibrator<5>& calibrator) noexcept
{
    // Right wing: b*(1 + rho) <= 2
    calibrator.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) -> double
        {
            // x = [a, b, rho, m, sigma]
            const double b{ x[1] };
            const double rho{ x[2] };

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
        nullptr
    );

    // Left wing: b*(1 - rho) <= 2
    calibrator.addInequalityConstraint(
        +[](unsigned, const double* x, double* grad, void*) -> double
        {
            const double b{ x[1] };
            const double rho{ x[2] };

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
        nullptr
    );
}

void SVI::addCalendarConstraint(Calibrator<5>& calibrator,
    std::vector<ConstraintCtx>& contexts) noexcept
{
    // For each k, add one no calendar arbitrage constraint
    for (auto& ctx : contexts)
    {
        calibrator.addInequalityConstraint(
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
                    grad[3] = b * (rho + xi * invR);   // -∂w/∂m
                    grad[4] = -b * (sigma * invR);      // -∂w/∂σ
                }

                // NLopt expects c(x) ≤ 0.
                // Here:  c(x) = prevWk + eps − w(k)
                // Enforces w(k) ≥ prevWk + eps  (no calendar arbitrage, with small tolerance)
                return c.prevWk + c.eps - wK;
            },
            &ctx
        );
    }
}

void SVI::addConvexityConstraint(Calibrator<5>& calibrator, std::vector<double>& kStorage) noexcept
{
    // For each k, add one convexity constraint g(k) ≥ 0
    for (std::size_t i = 0; i < kStorage.size(); ++i)
    {
        calibrator.addInequalityConstraint(
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
                const double g{ SVI::gk(a, b, rho, m, sigma, k, p) };

                if (grad)
                {
                    // ∇g = [dg/da, dg/db, dg/dρ, dg/dm, dg/dσ]
                    const auto dg{ SVI::gkGrad(a, b, rho, m, sigma, k, p) };
                    for (unsigned j = 0; j < 5 && j < n; ++j)
                        grad[j] = -dg[j]; // c(x) = -g(k)
                }

                // NLopt enforces c(x) ≤ tol. 
                // Here c(x) = -g(k), so -g(k) ≤ tol  
                // Enforces g(k) ≥ -tol  (≈ g(k) ≥ 0 with small slack)
                return -g;
            },
            &kStorage[i]  // pointer to stable local copy
        );
    }
}

void SVI::setMinObjective(Calibrator<5>& calibrator, const ObjCtx& obj) noexcept
{
    calibrator.setMinObjective(
        +[](unsigned n, const double* x, double* grad, void* data) -> double
        {
            // Reinterpret opaque pointer
            const ObjCtx& obj{ *static_cast<const ObjCtx*>(data) };

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

            for (std::size_t i = 0; i < obj.n; ++i)
            {
                // Extract variables
                const double k{ obj.k[i] };
                const double wKMarket{ obj.wK[i] };

                // Precompute variables
                const double xi{ k - m };
                const double R{ std::hypot(xi, sigma) };
                const double wKModel{ std::fma(b, rho * xi + R, a) };

                // Calculate SSE between model and market variance
                const double r{ wKModel - wKMarket };
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
        const_cast<ObjCtx*>(&obj) // nlopt needs void*
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

std::array<double, 5> SVI::gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
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

void SVI::evalCal(const SVISlice& sviSlice,
    const Calibrator<5>& calibrator,
    const std::vector<double>& kSlice,
    const std::vector<double>& wKPrevSlice) noexcept
{
    // Extract parameters
    const double a{ sviSlice.a };
    const double b{ sviSlice.b };
    const double rho{ sviSlice.rho };
    const double m{ sviSlice.m };
    const double sigma{ sviSlice.sigma };
 
    // Extract config attributes
    const double eps{ calibrator.eps()};
    const double tol{ calibrator.tol()};

    // ---- wMin constraint violation ----
    const double S = std::sqrt(std::max(0.0, 1.0 - rho * rho));
    const double wMin = std::fma(b, sigma * S, a);

    UV_WARN(eps - wMin > tol,
        std::format("No-arbitrage: wMin violated: wMin = {:.6f} < eps = {:.6f}",
            wMin, eps));

    // ---- Lee wing-slope constraint violations ----
    const double leeRight = b * (1.0 + rho) - 2.0;
    const double leeLeft = b * (1.0 - rho) - 2.0;

    UV_WARN(leeRight > tol,
        std::format("No-arbitrage: right wing slope > 2: b*(1+rho) = {:.6f} "
            "(b = {:.6f}, rho = {:.6f})",
            b * (1.0 + rho), b, rho));

    UV_WARN(leeLeft > tol,
        std::format("No-arbitrage: left wing slope > 2: b*(1-rho) = {:.6f} "
            "(b = {:.6f}, rho = {:.6f})",
            b * (1.0 - rho), b, rho));

    // ---- Convexity g(k) violation ----
    double gMin = std::numeric_limits<double>::infinity();
    double kAtMin = 0.0;
    std::size_t nViol = 0;

    for (double k : kSlice)
    {
        const GKPrecomp pre{ a, b, rho, m, sigma, k };
        const double g = SVI::gk(a, b, rho, m, sigma, k, pre);
        if (g < gMin) { gMin = g; kAtMin = k; }
        if (g < -tol) ++nViol;
    }

    UV_WARN(gMin < -tol,
        std::format("No-arbitrage: convexity violated: min g(k) = {:.6e} at k = {:.4f} "
            "(violations {} / {}, tol = {:.2e})",
            gMin, kAtMin, nViol, kSlice.size(), tol));

    // ---- Calendar no-arb: w_curr(k) >= w_prev(k) + eps ----
    std::size_t nCalViol = 0;
    double minMargin = std::numeric_limits<double>::infinity();
    double kAtWorst = 0.0;

    for (std::size_t i = 0; i < kSlice.size(); ++i)
    {
        const double k = kSlice[i];
        const double wPrev = wKPrevSlice[i] + eps;
        const double wCurr = SVI::wk(a, b, rho, m, sigma, k);
        const double margin = wCurr - wPrev; // should be >= 0
        if (margin < minMargin) { minMargin = margin; kAtWorst = k; }
        if (margin < -tol) ++nCalViol;
    }

    UV_WARN(minMargin < -tol,
        std::format("No-arbitrage: calendar violated: "
            "min[w_curr - (w_prev+eps)] = {:.6e} at k = {:.4f} "
            "(violations {} / {}, tol = {:.2e}, eps = {:.2e})",
            minMargin, kAtWorst, nCalViol, kSlice.size(),
            tol, eps));
}
  
std::vector<double> SVI::makewKSlice(const std::vector<double>& kSlice,
    double a, double b, double rho, double m, double sigma) noexcept
{
    std::vector<double> wKSlice;
    wKSlice.reserve(kSlice.size());

    std::transform(
        kSlice.begin(), kSlice.end(),
        std::back_inserter(wKSlice),
        [a, b, rho, m, sigma](double k) noexcept 
        {
            return wk(a, b, rho, m, sigma, k);
        });

    return wKSlice;
}