// SPDX-License-Identifier: Apache-2.0

#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(IntegrationSVISurfaceInvariants, RawSlicesHavePositiveVarianceAndButterflyDiagnostic)
{
    const uv::Vector<uv::models::svi::Params<double>> params{
        {0.5, 0.030, 0.20, -0.35, 0.02, 0.30},
        {1.0, 0.055, 0.24, -0.30, 0.02, 0.34},
        {2.0, 0.100, 0.28, -0.25, 0.02, 0.38}
    };
    const std::vector<double> logKF{-0.40, -0.25, -0.10, 0.0, 0.10, 0.25, 0.40};

    for (const auto& p : params)
    {
        for (double k : logKF)
        {
            const double w =
                uv::models::svi::totalVariance(p.a, p.b, p.rho, p.m, p.sigma, k);
            const double g = uv::models::svi::gk(p.a, p.b, p.rho, p.m, p.sigma, k);

            EXPECT_TRUE(std::isfinite(w)) << "t=" << p.t << " k=" << k;
            EXPECT_TRUE(std::isfinite(g)) << "t=" << p.t << " k=" << k;
            EXPECT_GT(w, 0.0) << "t=" << p.t << " k=" << k;
            EXPECT_GT(g, 0.0) << "t=" << p.t << " k=" << k;
        }
    }
}

TEST(IntegrationSVISurfaceInvariants, TotalVarianceIsCalendarMonotoneOnSharedLogMoneyness)
{
    const uv::Vector<uv::models::svi::Params<double>> params{
        {0.5, 0.030, 0.20, -0.35, 0.02, 0.30},
        {1.0, 0.055, 0.24, -0.30, 0.02, 0.34},
        {2.0, 0.100, 0.28, -0.25, 0.02, 0.38}
    };
    const std::vector<double> logKF{-0.40, -0.25, -0.10, 0.0, 0.10, 0.25, 0.40};

    for (double k : logKF)
    {
        double previous = 0.0;
        for (const auto& p : params)
        {
            const double w =
                uv::models::svi::totalVariance(p.a, p.b, p.rho, p.m, p.sigma, k);
            EXPECT_GE(w + 1e-14, previous) << "k=" << k << " t=" << p.t;
            previous = w;
        }
    }
}

TEST(IntegrationSVISurfaceInvariants, BuiltSurfaceHasFinitePositiveVols)
{
    const std::vector<double> maturities{0.5, 1.0, 2.0};
    const std::vector<double> forwards{100.0, 101.0, 102.0};
    const std::vector<double> strikes{70.0, 85.0, 100.0, 115.0, 130.0};
    const std::vector<double> moneyness{0.70, 0.85, 1.00, 1.15, 1.30};
    const uv::core::Matrix<double> inputVol{maturities.size(), strikes.size(), 0.25};
    const uv::core::VolSurface<double>
        input{maturities, forwards, strikes, moneyness, inputVol};
    const uv::Vector<uv::models::svi::Params<double>> params{
        {0.5, 0.030, 0.20, -0.35, 0.02, 0.30},
        {1.0, 0.055, 0.24, -0.30, 0.02, 0.34},
        {2.0, 0.100, 0.28, -0.25, 0.02, 0.38}
    };

    const auto surface = uv::models::svi::buildSurface(input, params);

    ASSERT_EQ(surface.numMaturities(), maturities.size());
    ASSERT_EQ(surface.numStrikes(), strikes.size());
    for (std::size_t i = 0; i < surface.numMaturities(); ++i)
    {
        for (std::size_t j = 0; j < surface.numStrikes(); ++j)
        {
            EXPECT_TRUE(std::isfinite(surface.vol()[i][j]));
            EXPECT_GT(surface.vol()[i][j], 0.0);
        }
    }
}
