/**
* GaussLegendre.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Quadrature/TanHSinH.hpp"
#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"
#include <boost/math/special_functions/lambert_w.hpp>
#include <utility>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <numbers>

namespace uv
{
    TanHSinH::TanHSinH(const unsigned int N) :
        N_(static_cast<size_t>(N)),
        h_
        (
            // Define optimal step size using a heuristic rule
            boost::math::lambert_w0
            (
                2.0L * std::numbers::pi_v<long double> *static_cast<long double>(N)
            )
            / static_cast<long double>(N)
        )
    {
        // Resize M vector
        nodes_.resize(N_ + 1);

        for (unsigned int n = 0; n <= N_; ++n)
        {
            // Calculate n*h
            const long double nh{ static_cast<long double>(n) * h_ };

            // Calculate and store node values
            nodes_[n] = generateNode(nh);
        }
    }

    TanHSinH::Node TanHSinH::generateNode(const long double nh) noexcept
    {
        // Calculate qn term
        long double qn{ std::exp(-std::numbers::pi_v<long double> *std::sinh(nh)) };

        // Precompute repetitive calculations
        long double qnInv{ 1.0L / (1.0L + qn) };

        // Calculate yn term
        long double yn{ (2.0L * qn * qnInv) };

        // Calculate and return Node struct
        return
        {
            yn,                                                                 // yn term
            1.0L - yn,                                                          // Abscissas value
            qnInv * yn * std::numbers::pi_v<long double> *std::cosh(nh)        // Weight value
        };
    }

    void TanHSinH::printGrid() const noexcept
    {
        constexpr int idx_w = 6;
        constexpr int col_w = 24;

        std::ostringstream oss;

        // Left align the title
        oss << std::left << "\nFixed Tanh-Sinh Grid\n";

        // Reset to default (right) for numeric columns
        oss << std::right;

        // Header
        oss << std::setw(col_w) << "x_n (node)" << ' '
            << std::setw(col_w) << "w_n (weight)" << '\n';

        // Separator
        oss << std::string(idx_w + 1 + col_w + 1 + col_w, '-') << '\n';

        // Body
        oss << std::scientific << std::setprecision(16);
        for (std::size_t n = 0; n < N_; ++n)
        {
            oss << std::setw(idx_w) << std::left << n
                << std::setw(col_w) << std::right << nodes_[n].x << ' '
                << std::setw(col_w) << std::right << nodes_[n].w << '\n';
        }

        UV_INFO(oss.str());
    }
}