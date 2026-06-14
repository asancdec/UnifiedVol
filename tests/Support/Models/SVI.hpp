// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Matrix.hpp"
#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"
#include "Support/Diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace uv::tests::models::svi
{
struct SyntheticSurfaceCase
{
    std::vector<double> maturities;
    std::vector<double> logKF;
    Vector<uv::models::svi::Params<double>> truth;
    core::Matrix<double> logKFMatrix;
    core::Matrix<double> totalVariance;
};

inline SyntheticSurfaceCase makeSyntheticSurfaceCase(
    std::vector<double> maturities,
    std::vector<double> logKF,
    Vector<uv::models::svi::Params<double>> truth
)
{
    const std::size_t rows{maturities.size()};
    const std::size_t cols{logKF.size()};
    SyntheticSurfaceCase data{
        .maturities = std::move(maturities),
        .logKF = std::move(logKF),
        .truth = std::move(truth),
        .logKFMatrix = core::Matrix<double>{rows, cols},
        .totalVariance = core::Matrix<double>{rows, cols}
    };

    for (std::size_t i = 0; i < data.maturities.size(); ++i)
    {
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            data.logKFMatrix[i][j] = data.logKF[j];
            const auto& p = data.truth[i];
            data.totalVariance[i][j] = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
        }
    }

    return data;
}

inline SyntheticSurfaceCase makeDefaultSyntheticSurfaceCase()
{
    return makeSyntheticSurfaceCase(
        {0.5, 1.0},
        {-0.25, -0.10, 0.0, 0.10, 0.25},
        {{0.5, 0.025, 0.30, -0.35, 0.02, 0.30}, {1.0, 0.040, 0.35, -0.25, 0.04, 0.35}}
    );
}

inline SyntheticSurfaceCase makePerformanceSyntheticSurfaceCase()
{
    return makeSyntheticSurfaceCase(
        {0.08, 0.25, 0.5, 1.0, 2.0},
        {-0.45, -0.30, -0.15, 0.0, 0.15, 0.30, 0.45},
        {{0.08, 0.010, 0.22, -0.70, -0.03, 0.22},
         {0.25, 0.025, 0.26, -0.65, -0.02, 0.26},
         {0.5, 0.045, 0.30, -0.55, -0.01, 0.30},
         {1.0, 0.080, 0.34, -0.45, 0.00, 0.36},
         {2.0, 0.140, 0.38, -0.35, 0.01, 0.44}}
    );
}

inline ErrorDiagnostics totalVarianceDiagnostics(
    const Vector<uv::models::svi::Params<double>>& calibrated,
    const SyntheticSurfaceCase& data
)
{
    double sumAbs{};
    double sumSquared{};
    double maxAbs{};
    std::size_t count{};

    for (std::size_t i = 0; i < calibrated.size(); ++i)
    {
        for (std::size_t j = 0; j < data.logKF.size(); ++j)
        {
            const auto& p = calibrated[i];
            const double fitted = uv::models::svi::totalVariance(
                p.a,
                p.b,
                p.rho,
                p.m,
                p.sigma,
                data.logKF[j]
            );
            const double error = std::abs(fitted - data.totalVariance[i][j]);

            sumAbs += error;
            sumSquared += error * error;
            maxAbs = std::max(maxAbs, error);
            ++count;
        }
    }

    return {
        sumAbs / static_cast<double>(count),
        maxAbs,
        std::sqrt(sumSquared / static_cast<double>(count))
    };
}
} // namespace uv::tests::models::svi
