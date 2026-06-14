// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Price/Pricer.hpp"
#include "Support/Tolerances.hpp"

#include <array>
#include <cmath>
#include <gtest/gtest.h>

namespace
{
template <typename F> double centralDifference(F&& f, double x, double h)
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}
} // namespace

TEST(RegressionHestonGradient, MatchesCentralFiniteDifference)
{
    uv::models::heston::price::Pricer<double, 160> pricer{};

    const double kappa = 2.0;
    const double theta = 0.05;
    const double sigma = 0.4;
    const double rho = -0.6;
    const double v0 = 0.04;
    const double t = 1.25;
    const double dF = 0.97;
    const double F = 100.0;
    const double K = 105.0;

    const auto analytic =
        pricer.callPriceWithGradient(kappa, theta, sigma, rho, v0, t, dF, F, K);

    const std::array<double, 5> params{kappa, theta, sigma, rho, v0};
    const std::array<double, 5> steps{1e-4, 1e-5, 1e-5, 1e-5, 1e-5};

    for (std::size_t idx = 0; idx < params.size(); ++idx)
    {
        const double finiteDifference = centralDifference(
            [&](double shifted)
            {
                auto p = params;
                p[idx] = shifted;
                return pricer.callPrice(p[0], p[1], p[2], p[3], p[4], t, dF, F, K);
            },
            params[idx],
            steps[idx]
        );

        EXPECT_NEAR(
            analytic[idx + 1],
            finiteDifference,
            uv::tests::tolerance::FiniteDifference
        ) << "parameter index "
          << idx;
    }
}
