// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Detail/Initialize.hpp"

namespace uv::models::svi::detail
{
std::array<double, 4> coldGuess() noexcept
{
    return {
        0.1,  // b: wing slope scale
        -0.5, // rho: skew correlation
        0.1,  // m: horizontal translation
        0.1   // sigma: curvature scale
    };
}

std::array<double, 4> warmGuess(const Params<double>& params) noexcept
{
    return {
        params.b,    // b
        params.rho,  // rho
        params.m,    // m
        params.sigma // sigma
    };
}

std::array<double, 4> lowerBounds(double logKFMin) noexcept
{
    return {
        0.001,           // b
        -0.9999,         // rho
        10.0 * logKFMin, // m
        0.01             // sigma
    };
}

std::array<double, 4> upperBounds(double logKFMax) noexcept
{
    return {
        2.0,             // b
        0.9999,          // rho
        10.0 * logKFMax, // m
        10.0             // sigma
    };
}

} // namespace uv::models::svi::detail
