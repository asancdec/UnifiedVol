/**
* SVI.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/

#include "VolModels/SVI.hpp"

#include <nlopt.hpp>

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


// Parameter calibration occurs in the constructor
SVI::SVI(const VolSurface& impVolSurf, unsigned long int maxEval) : modelVolSurf(impVolSurf)
{
    // Calculate log-forward moneyness (log(K/S))
    auto k = impVolSurf.k();

    // Calculate variance time implied surface
    auto wMkt = impVolSurf.wMkt();

    // Determine the ATM index and the indexes before and after
    auto [atmInd, prevInd, nextInd] = impVolSurf.findATMIndices();

    // Clear previous parameters
    this->rawParams.clear();

    // Calibrate parameters on a maturity basis
    for (size_t i = 0; i < impVolSurf.maturities.size(); ++i)
    {   
        double T{ impVolSurf.maturities[i] };

        // Initialize SliceData struct
        SliceData sliceData{ k[i], wMkt[i] };

        // Nelder–Mead / Subplex derivative-free optimization
        nlopt::opt opt(nlopt::LN_SBPLX, 5);

        // a     : total variance level (ATM)             -> must be non-negative, allows small variance
        // b     : slope of the wings                     -> must be positive to avoid sqrt of negative
        // rho   : correlation / skew parameter           -> constrained to (-1, 1) to prevent arbitrage
        // m     : horizontal shift of the smile          -> center of log-forward moneyness
        // sigma : curvature / width of the smile         -> must be positive to avoid sqrt of negative
        opt.set_lower_bounds({ 0.0, 1e-8,  -0.999,  -5.0, 1e-8 });
        opt.set_upper_bounds({ 5.0,  10.0,   0.999,   5.0, 5.0 });

        // Set minimization function
        opt.set_min_objective(sviObjective, &sliceData);

        // Set maximum number of evaluations
        opt.set_maxeval(maxEval);

        // Define initial guess
        std::vector<double> x = this->initGuess(impVolSurf, sliceData, T, i, atmInd, prevInd, nextInd);
        double minf;

        // Run optimization
        nlopt::result result{ opt.optimize(x, minf) };

        // Unpack the SVI parameters 
        double a = x[0], b = x[1], rho = x[2], m = x[3], sigma = x[4];

        // Store calibrated params
        this->rawParams[impVolSurf.maturities[i]] = { a, b, rho, m, sigma };
    }

    // Build model-implied volatilities using calibrated parameters 
    buildModelSurface(impVolSurf, k);
}

// Initial guess using SVI Jump-Wings (SVI-JW) parameterization
std::vector<double> SVI::initGuess(const VolSurface& impVolSurf, const SliceData& sliceData, const double T,
    size_t i, size_t atmInd, size_t prevInd, size_t nextInd) const
{
    // Convenient local variables
    const auto& k = sliceData.k;
    const auto& wMkt = sliceData.wMkt;
    const auto& vols = impVolSurf.vols[i];

    // Determine ATM total variance
    double w{ wMkt[atmInd] };

    // Calculate annualized ATM variance 
    double v{ wMkt[atmInd] / T };

    // Calculate ATM skew
    double psi{ (vols[nextInd] - vols[prevInd]) / (k[nextInd] - k[prevInd])};

    // Calculate slope of the left (put) wing
    double p{ (wMkt[atmInd] - wMkt[0]) / (k[atmInd] - k[0]) };

    // Calculate slope of the right (call) wing
    double c{ (wMkt.back() - wMkt[atmInd]) / (k.back() - k[atmInd]) };

    // Calculate minimum implied variance
    double vMin{ *std::min_element(wMkt.begin(), wMkt.end()) / T };

    // Calculate raw SVI parameters with the SVI JW estimates
    double b{ this->b(w, c, p) };
    double rho{ this->rho(p, w, b) };
    double beta{ std::clamp(this->beta(rho, psi, w, b), -0.9999, 0.9999)};
    double alpha{ this->alpha(beta) };
    double m{ this->m(v, vMin, T, b, rho, alpha) };
    double sigma{ this->sigma(m, alpha) };
    double a{ this->a(vMin, T, b, sigma, rho) };


    std::cout << "ATM total variance w: " << w << std::endl;
    std::cout << "Annualized ATM variance v: " << v << std::endl;
    std::cout << "ATM skew psi: " << psi << std::endl;
    std::cout << "Left wing slope p: " << p << std::endl;
    std::cout << "Right wing slope c: " << c << std::endl;
    std::cout << "Minimum implied variance vMin: " << vMin << std::endl;
    std::cout << "SVI parameter b: " << b << std::endl;
    std::cout << "SVI parameter rho: " << rho << std::endl;
    std::cout << "SVI parameter beta: " << beta << std::endl;
    std::cout << "SVI parameter alpha: " << alpha << std::endl;
    std::cout << "SVI parameter m: " << m << std::endl;
    std::cout << "SVI parameter sigma: " << sigma << std::endl;
    std::cout << "SVI parameter a: " << a << std::endl;


    return { a, b, rho, m, sigma };
}

double SVI::b(double w, double c, double p)
{
    return std::sqrt(w) / 2 * (c + p);
}

double SVI::rho(double p, double w, double b)
{
    return 1 - p * std::sqrt(w) / b;
}

double SVI::a(double vMin, double T, double b, double sigma, double rho)
{
    return vMin * T - b * sigma * std::sqrt(1 - rho * rho);
}

double SVI::m(double v, double vMin, double T, double b, double rho, double alpha)
{   
    int signAlpha{ std::signbit(alpha) ? -1 : 1 };
    return (v - vMin) * T / (b * (-rho + signAlpha * std::sqrt(1 + alpha * alpha) - alpha * std::sqrt(1 - rho * rho)));
}

double SVI::sigma(double m, double alpha)
{
    return alpha * m;
}

double SVI::beta(double rho, double psi, double w, double b)
{
    return rho - 2 * psi * std::sqrt(w) / b;
}

double SVI::alpha(double beta)
{
    int signBeta{ std::signbit(beta) ? -1 : 1 };
    return signBeta * std::sqrt(1 / (beta * beta) - 1);

}

// Build SVI volatility surface
void SVI::buildModelSurface(const VolSurface& impVolSurf, const std::vector<std::vector<double>>& k)
{      
    // Iterate by rows (maturities)
    for (size_t i = 0; i < impVolSurf.maturities.size(); ++i) 
    {
        double T{ impVolSurf.maturities[i] };
        const auto& params = this->rawParams[T];

        // Calculate SVI volatility at each maturity and strike
        for (size_t j = 0; j < k[i].size(); ++j) 
        {
            double wModel{ sviVariance(params.a, params.b, params.rho, params.m, params.sigma, k[i][j]) };
            this->modelVolSurf.vols[i][j] = std::sqrt(wModel / T);
        }
    }
}

// Static objective function
double SVI::sviObjective(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    // Cast the void pointer to SliceData (market slice for one maturity)
    auto* slice{ reinterpret_cast<SliceData*>(data) };
    double error{ 0.0 };

    // Loop over strikes / log-forward moneyness
    for (size_t i = 0; i < slice->k.size(); ++i)
    {
        double k{ slice->k[i] };      // current log-forward moneyness
        double wMkt{ slice->wMkt[i] }; // market total variance at this strike

        // Unpack the SVI parameters 
        double a = x[0], b = x[1], rho = x[2], m = x[3], sigma = x[4];

        // Compute model total variance at this strike
        double wModel{ SVI::sviVariance(a, b, rho, m, sigma, k) };

        // Squared error contribution
        double diff{ wModel - wMkt };
        error += diff * diff;
    }

    return error; // return total sum of squared errors for the slice
}

// Static model variance calculation function
double SVI::sviVariance(double a, double b, double rho, double m, double sigma, double k)
{
    return a + b * (rho * (k - m) + std::sqrt((k - m) * (k - m) + sigma * sigma));
}

