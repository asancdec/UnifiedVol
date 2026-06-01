// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Detail/Initialize.hpp"

namespace uv::models::svi::detail
{
std::array<double, 4> coldGuess() noexcept
{
    return {0.1, -0.5, 0.1, 0.1};
}

std::array<double, 4> warmGuess(const Params<double>& params) noexcept
{
    return {

        params.b,
        params.rho,
        params.m,
        params.sigma
    };
}

std::array<double, 4> lowerBounds(double logKFMin) noexcept
{
    return {0.001, -0.9999, 10.0 * logKFMin, 0.01};
}

std::array<double, 4> upperBounds(double logKFMax) noexcept
{
    return {2.0, 0.9999, 10.0 * logKFMax, 10.0};
}

} // namespace uv::models::svi::detail