// SPDX-License-Identifier: Apache-2.0

#include <array>

namespace uv::models::heston::calibrate::detail
{

[[nodiscard]] constexpr std::array<double, 5> initGuess() noexcept
{
    return {
        2.5,   // kappa: mean-reversion speed
        0.09,  // theta: long-run variance
        1.0,   // sigma: volatility of variance
        -0.75, // rho: spot/variance correlation
        0.20   // v0: initial variance
    };
}

[[nodiscard]] constexpr std::array<double, 5> lowerBounds() noexcept
{
    return {
        0.001,  // kappa
        0.001,  // theta
        0.001,  // sigma
        -0.999, // rho
        0.001   // v0
    };
}

[[nodiscard]] constexpr std::array<double, 5> upperBounds() noexcept
{
    return {
        15.0,  // kappa
        0.5,   // theta
        10.0,  // sigma
        0.999, // rho
        0.5    // v0
    };
}

template <typename Policy> void setGuessBounds(opt::ceres::Optimizer<Policy>& optimizer)
{
    optimizer.initialize(initGuess(), lowerBounds(), upperBounds());
}

} // namespace uv::models::heston::calibrate::detail
