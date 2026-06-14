// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Price/Pricer.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

TEST(Heston, IntrinsicValue)
{
    using Real = long double;

    uv::models::heston::price::Pricer<Real> hestonPricer{};

    Real t = 1e-14;
    Real dF = 0.98;
    Real F = 1.0;
    Real K = 0.98;

    Real kappa = 4.0;
    Real theta = 0.052;
    Real sigma = 1e-14;
    Real rho = -0.7;
    Real v0 = 0.22;

    Real callHeston = hestonPricer.callPrice(kappa, theta, sigma, rho, v0, t, dF, F, K);

    Real intrinsic = dF * std::max(F - K, Real{0});

    EXPECT_NEAR(callHeston, intrinsic, 1e-11);
}
