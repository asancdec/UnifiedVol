// SPDX-License-Identifier: Apache-2.0

#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Params.hpp"

#include <gtest/gtest.h>
#include <vector>

TEST(IntegrationModelsSVIBuildSurface, BuildsVolSurfaceFromKnownParams)
{
    const std::vector<double> maturities{1.0, 2.0};
    const std::vector<double> forwards{100.0, 100.0};
    const std::vector<double> strikes{90.0, 100.0, 110.0};
    const std::vector<double> moneyness{0.9, 1.0, 1.1};
    const uv::core::Matrix<double> inputVol{2, 3, 0.2};
    const uv::core::VolSurface<double>
        input{maturities, forwards, strikes, moneyness, inputVol};
    const uv::Vector<uv::models::svi::Params<double>> params{
        {1.0, 0.02, 0.10, 0.0, 0.0, 0.20},
        {2.0, 0.02, 0.10, 0.0, 0.0, 0.20}
    };

    const auto output = uv::models::svi::buildSurface(input, params);

    EXPECT_EQ(output.numMaturities(), input.numMaturities());
    EXPECT_EQ(output.numStrikes(), input.numStrikes());
    EXPECT_DOUBLE_EQ(output.maturities()[1], input.maturities()[1]);
    EXPECT_DOUBLE_EQ(output.strikes()[2], input.strikes()[2]);
    EXPECT_GT(output.vol()[0][0], 0.0);
    EXPECT_GT(output.vol()[1][2], 0.0);
}
