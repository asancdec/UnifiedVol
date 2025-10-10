/**
* StandardGL.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Quadrature/StandardGL.hpp"
#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"

#include <Eigen/Dense>
#include <utility>
#include <iomanip>
#include <sstream>


using uv::ErrorCode;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::VectorXd;

StandardGL::StandardGL(const int N, const double alpha) : N_(N)
{

    // Input validation
    UV_REQUIRE(N_ > 0, ErrorCode::InvalidArgument, "Gauss-Laguerre: N must be > 0");
    UV_REQUIRE(alpha > -1.0, ErrorCode::InvalidArgument, "Gauss-Laguerre: alpha must be > -1");

    // Build Jacobi (symmetric tridiagonal) matrix J
    // Diagonal:   a_k = 2k + alpha + 1,   k = 0..N-1
    MatrixXd J = MatrixXd::Zero(N_, N_);
    for (int k = 0; k < N_; ++k)
    {
        J(k, k) = 2.0 * k + alpha + 1.0;
    }

    // Off-diag:   b_k = sqrt(k*(k+alpha)), k = 1..N-1
    for (int k = 1; k < N_; ++k)
    {
        double bk{ std::sqrt(k * (k + alpha)) };
        J(k - 1, k) = bk;
        J(k, k - 1) = bk;
    }

    // Eigen-decomposition (nodes = eigenvalues, weights from eigenvectors) 
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(J);
    if (es.info() != Eigen::Success) 
    {
        uv::raise(ErrorCode::LinearAlgebra,
            "Gauss-Laguerre: eigen decomposition failed (SelfAdjointEigenSolver)");
    }

    const VectorXd evals{ es.eigenvalues() };     // x_i (ascending)
    const MatrixXd V{ es.eigenvectors() };        // columns are normalized eigenvectors v_i

    // Mu is defined as the zeroth moment of the weight function
    // For Gauss Laguerre:
    // w(x) = x^alpha * exp(-x)
    // x E [0,inf)
    const double mu0{ std::tgamma(alpha + 1.0) };  // mu0 = ∫_0^inf w(x) dx = Gamma(alpha+1)

    xk_.resize(N_);
    wk_.resize(N_);

    for (int i = 0; i < N_; ++i)
    {
        xk_[i] = evals(i);         // Quadrature node (root of Laguerre polynomial)
        double v0{ V(0, i) };      // First component of eigenvector i
        wk_[i] = mu0 * (v0 * v0);  // Quadrature weight corresponding to xk_[i]
    }
}

void StandardGL::printGrid() const noexcept
{
    std::ostringstream oss;
    oss << "\n";
    oss << std::fixed << std::setprecision(12);

    oss << "Gauss-Laguerre grid (N=" << N_ << ")\n";
    oss << std::left << std::setw(6) << "#"
        << std::right << std::setw(22) << "x_i (node)"
        << std::setw(22) << "w_i (weight)" << '\n';
    oss << std::string(6 + 22 + 22, '-') << '\n';

    for (int i = 0; i < N_; ++i)
    {
        oss << std::left << std::setw(6) << i
            << std::right << std::scientific << std::setprecision(12)
            << std::setw(22) << xk_[i]
            << std::setw(22) << wk_[i]
            << '\n';
    }
    UV_INFO(oss.str());
}