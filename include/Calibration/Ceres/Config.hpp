/**
* Config.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <array>    
#include <string_view>
#include <cstddef>

namespace uv::cal::ceres
{
	template <std::size_t N>
	struct Config
    {
        unsigned maxEval;                               // Maximum number of iterations / function evaluations
        double functionTol;                             // Function tolerance → stop when relative cost improvement < threshold
        double paramTol;                                // Parameter tolerance → stop when parameter updates are below threshold
        std::array<std::string_view, N> paramNames;     // Parameter names (for logging and diagnostics)
        double gradientTol{ double(0.0) };              // Gradient tolerance → stop when ∥Jᵀr∥ < threshold (stationarity)
        double lossScale{ double(1.0) };                // Loss function scale (e.g., δ for Huber or Cauchy)
        bool verbose{false};                            // Logs the full Ceres calibration report 
    };
}
