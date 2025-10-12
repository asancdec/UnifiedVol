/**
* GaussLegendre.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Quadrature/GaussLaguerre.hpp"
#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"

#include <Eigen/Dense>
#include <utility>
#include <iomanip>
#include <sstream>
#include <cmath>

using uv::ErrorCode;

GaussLaguerre::GaussLaguerre(const int N, const double alpha) : N_(N), alpha_(static_cast<long double>(alpha))
{

    // Input validation
    UV_REQUIRE(N_ > 0, ErrorCode::InvalidArgument, "GL: N must be > 0");
    UV_REQUIRE(alpha_ > -1.0, ErrorCode::InvalidArgument, "GL: alpha must be > -1");

    // Jacobi matrix coefficients: diagonal a_i, subdiagonal b_i
    Eigen::VectorXd d(N_), e(N_ - 1);   
    for (int i = 0; i < N_; ++i) d(i) = 2.0 * i + static_cast<double>(alpha_) + 1.0;
    for (int i = 1; i < N_; ++i) e(i - 1) = std::sqrt(i * (i + static_cast<double>(alpha_)));

    // Eigen-decomposition (nodes = eigenvalues)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.computeFromTridiagonal(d, e);

    UV_REQUIRE(es.info() == Eigen::Success, ErrorCode::LinearAlgebra, "GL: eigen decomposition failed" );

    // Extract nodes and compute weights
    const Eigen::VectorXd& evals{ es.eigenvalues() };
    xk_.resize(N_);
    wk_.resize(N_);
    for (int i = 0; i < N_; ++i)
    {   
        const long double x{ static_cast<long double>(evals(i)) };
        xk_[i] = x;        
        wk_[i] = weight(x);
    }
}

long double GaussLaguerre::Lprime(const long double xi) const noexcept
{
    // Analytic identity:
    // d/dx L_N^(α)(x) = -L_{N-1}^(α+1)(x)

    if (N_ == 0) return 0.0L;   // d/dx L_0 = 0
    if (N_ == 1) return -1.0L;  // d/dx(1 + α − x) = −1

    const long double a{ alpha_ + 1.0L };  // shifted α
    const int n{ N_ - 1 };                 // degree for the inner Laguerre

    // Recurrence for L_{n}^{(a)}(x)
    long double Lnm2{ 1.0L };               // L_0^(a)
    long double Lnm1{ 1.0L + a - xi };      // L_1^(a)
    long double Ln{ 0.0L };

    for (int k = 2; k <= n; ++k)
    {
        const long double kk{ static_cast<long double>(k) };
        Ln = ((2.0L * kk - 1.0L + a - xi) * Lnm1 - (kk - 1.0L + a) * Lnm2) / kk;
        Lnm2 = Lnm1;
        Lnm1 = Ln;
    }
    return -Lnm1; // L'_N^(α)(x) = −L_{N−1}^(α+1)(x)
}

long double GaussLaguerre::weight(const long double xi) const noexcept
{   
    // Formula:
     //    w_i = Γ(N+α+1) / [ Γ(N+1) · x_i · (L'_N^(α)(x_i))² ]

    const long double dL{ Lprime(xi) };        // Derivative of L_N^(α) at x_i
    if (dL == 0.0L) return 0.0L;              // safeguard 

    // Trasnform to log-domain to improve numerical stability
    // logΓ(N+α+1) − logΓ(N+1) − log(x_i) − 2·log|L'_N^(α)(x_i)|
    const long double logNum{ std::lgamma(static_cast<long double>(N_) + alpha_ + 1.0L) };
    const long double logDen{ std::lgamma(static_cast<long double>(N_) + 1.0L)
        + std::log(xi)
        + 2.0L * std::log(std::fabs(dL)) };

    return std::exp(logNum - logDen);
}

void GaussLaguerre::printGrid() const noexcept 
{
    constexpr int idx_w = 6, col_w = 24;
    std::ostringstream oss;
    oss << "\nGauss-Laguerre grid (N=" << N_ << ")\n";
    oss << std::left << std::setw(idx_w) << "#"
        << std::right << std::setw(col_w) << "x_i (node)" << ' '
        << std::right << std::setw(col_w) << "w_i (weight)" << '\n';
    oss << std::string(idx_w + 1 + col_w + 1 + col_w, '-') << '\n';

    oss << std::scientific << std::setprecision(16);
    for (int i = 0; i < N_; ++i) {
        oss << std::left << std::setw(idx_w) << i
            << std::right << std::setw(col_w) << xk_[i] << ' '
            << std::right << std::setw(col_w) << wk_[i] << '\n';
    }
    UV_INFO(oss.str());
}