/**
* SVI.cpp
* Author: Alvaro Sanchez de Carlos
*/


#include "Models/SVI/SVI.hpp"   

#include <algorithm>   
#include <cmath>      
#include <format>     
#include <iterator>   
#include <vector>

namespace uv
{
    ::std::array<double, 5> SVI::initGuess(const SliceData& slice) noexcept
    {
        const double b{ 0.1 };
        const double rho{ -0.5 };
        const double m{ 0.1 };
        const double sigma{ 0.1 };
        const double a{ slice.atmWT() - b * (-rho * m + ::std::hypot(m, sigma)) };

        return { a, b, rho, m, sigma };
    }

    ::std::array<double, 5> SVI::lowerBounds(const SliceData& slice) noexcept
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

    ::std::array<double, 5> SVI::upperBounds(const SliceData& slice) noexcept
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

    double SVI::wk(double a, double b, double rho, double m, double sigma, double k) noexcept
    {
        const double x{ k - m };                                     // x := k-m
        return ::std::fma(b, (rho * x + ::std::hypot(x, sigma)), a);    // w(k) = a + b*(rho*x + sqrt(x^2 + sigma^2)) 
    }

    double SVI::gk(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
    {
        return (p.A * p.A) - 0.25 * p.wkD1Squared * p.B + p.wkD2 / 2.0;  // g(k) = A^2 - B * (w')^2 / 4 + w''/2
    }

    ::std::array<double, 5> SVI::gkGrad(double a, double b, double rho, double m, double sigma, double k, const GKPrecomp& p) noexcept
    {
        // Precompute variables
        double invR5{ p.invRCubed * p.invR * p.invR };     // 1/R^5
        double wkSquaredInv{ 1.0 / (p.wk * p.wk) };        // 1/w^2

        // Calculate partial derivatives of g(k) with respect to w, w', and w''
        double dgdw{ p.A * k * p.wkD1 * wkSquaredInv + 0.25 * p.wkD1Squared * wkSquaredInv };  // ∂g/∂w   = A*k*w'/w^2 + (w')^2 / (4 * w^2)
        double dgdw1{ -p.A * k / p.wk - 0.5 * p.wkD1 * p.B };                                  // ∂g/∂w'  = -A*k/w - w'*B/2
        double dgdw2{ 0.5 };                                                                   // ∂g/∂w'' = 1/2

        // ---------- ∂w/∂θ ---------
        ::std::array<double, 5> dw
        {
            1.0,                         // ∂w/∂a   = 1
            rho * p.x + p.R,             // ∂w/∂b   = ρ (k-m) + R
            b * p.x,                     // ∂w/∂ρ   = b (k-m)
            -b * (rho + p.x * p.invR),   // ∂w/∂m   = -b ( ρ + (k-m)/R )
            b * sigma * p.invR           // ∂w/∂σ   = b (σ / R)
        };

        // ---------- ∂w′/∂θ ------
        ::std::array<double, 5> dw1
        {
            0.0,                                // ∂w′/∂a   = 0
            rho + p.x * p.invR,                 // ∂w′/∂b   = ρ + (k-m)/R
            b,                                  // ∂w′/∂ρ   = b
            -b * p.sigmaSquared * p.invRCubed,  // ∂w′/∂m   = -b * σ^2 / R^3
            -b * p.x * sigma * p.invRCubed      // ∂w′/∂σ   = -b * (k-m) * σ / R^3
        };

        // ---------- ∂w″/∂θ ----------
        ::std::array<double, 5> dw2
        {
            0.0,                                                                 // ∂w″/∂a   = 0
            p.sigmaSquared * p.invRCubed,                                        // ∂w″/∂b   = σ^2 / R^3
            0.0,                                                                 // ∂w″/∂ρ   = 0
            3.0 * b * p.sigmaSquared * p.x * invR5,                              // ∂w″/∂m   = 3 b σ^2 (k-m) / R^5
            b * (2 * sigma * p.invRCubed - 3.0 * p.sigmaSquared * sigma * invR5) // ∂w″/∂σ   = b( 2σ/R^3 - 3σ^3/R^5 )
        };

        // Chain rule: ∇g = (∂g/∂w)∇w + (∂g/∂w1)∇w1 + (∂g/∂w2)∇w2
        ::std::array<double, 5> dg{};
        for (int j = 0; j < 5; ++j)
        {
            dg[j] = dgdw * dw[j] + dgdw1 * dw1[j] + dgdw2 * dw2[j];
        }
        return dg;
    }

    ::std::vector<double> SVI::makewKSlice(const ::std::vector<double>& kSlice,
        double a, double b, double rho, double m, double sigma) noexcept
    {
        ::std::vector<double> wKSlice;
        wKSlice.reserve(kSlice.size());

        ::std::transform(
            kSlice.begin(), kSlice.end(),
            ::std::back_inserter(wKSlice),
            [a, b, rho, m, sigma](double k) noexcept
            {
                return wk(a, b, rho, m, sigma, k);
            });

        return wKSlice;
    }
}