// SPDX-License-Identifier: Apache-2.0

#include "Math/Functions/Black.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <array>
#include <cmath>
#include <gtest/gtest.h>

TEST(Heston, SigmaZeroMatchesBlack)
{
    using Real = long double;

    uv::models::heston::price::Pricer<Real> hestonPricer{};

    Real t = 3.0;
    Real dF = 0.98;
    Real F = 1.0;
    Real K = 0.98;

    Real kappa = 4.0;
    Real theta = 0.052;
    Real sigma = 1e-14;
    Real rho = -0.5;
    Real v0 = 0.22;

    Real callHeston = hestonPricer.callPrice(kappa, theta, sigma, rho, v0, t, dF, F, K);

    Real I = theta * t + (v0 - theta) * (1.0 - std::exp(-kappa * t)) / kappa;

    Real volBS = std::sqrt(I / t);

    Real callBS = uv::math::black::priceB76(t, dF, F, volBS, K);

    EXPECT_NEAR(callHeston, callBS, 1e-14);
}

TEST(RegressionHestonBlackLimit, MatchesBlackAcrossGridAndIsRhoIndependent)
{
    using Real = long double;

    uv::models::heston::price::Pricer<Real, 192> hestonPricer{};

    constexpr Real kappa = 4.0L;
    constexpr Real theta = 0.052L;
    constexpr Real sigma = 1e-14L;
    constexpr Real v0 = 0.22L;
    constexpr Real dF = 0.98L;
    constexpr Real F = 1.0L;

    const std::array<Real, 3> maturities{0.25L, 1.0L, 3.0L};
    const std::array<Real, 5> strikes{0.80L, 0.95L, 1.0L, 1.05L, 1.20L};
    const std::array<Real, 3> rhos{-0.90L, 0.0L, 0.90L};

    for (Real t : maturities)
    {
        const Real integratedVariance =
            theta * t + (v0 - theta) * (1.0L - std::exp(-kappa * t)) / kappa;
        const Real blackVol = std::sqrt(integratedVariance / t);

        for (Real K : strikes)
        {
            const Real black = uv::math::black::priceB76(t, dF, F, blackVol, K);
            const Real referenceHeston =
                hestonPricer
                    .callPrice(kappa, theta, sigma, rhos.front(), v0, t, dF, F, K);

            EXPECT_NEAR(referenceHeston, black, 1e-11)
                << "t=" << static_cast<double>(t) << " K=" << static_cast<double>(K);

            for (Real rho : rhos)
            {
                const Real heston =
                    hestonPricer.callPrice(kappa, theta, sigma, rho, v0, t, dF, F, K);

                EXPECT_NEAR(heston, referenceHeston, 1e-12)
                    << "t=" << static_cast<double>(t) << " K=" << static_cast<double>(K)
                    << " rho=" << static_cast<double>(rho);
            }
        }
    }
}
