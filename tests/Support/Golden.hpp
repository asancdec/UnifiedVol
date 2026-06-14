// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/SVI/Params.hpp"

#include <cstddef>
#include <filesystem>

namespace uv::tests::golden
{
struct VolPoint
{
    std::size_t maturity{};
    std::size_t strike{};
    double value{};
};

struct ExamplePipeline
{
    double tolerance{};
    Vector<models::svi::Params<double>> sviParams;
    models::heston::Params<double> hestonParams{0.0, 0.0, 0.0, 0.0, 0.0};
    Vector<VolPoint> marketVols;
    Vector<VolPoint> sviVols;
    Vector<VolPoint> hestonVols;
    double meanAbsVolError{};
    double maxAbsVolError{};
    double rmseVolError{};
};

struct SyntheticSVICalibration
{
    double tolerance{};
    Vector<models::svi::Params<double>> calibratedParams;
    double meanAbsTotalVarianceError{};
    double maxAbsTotalVarianceError{};
    double rmseTotalVarianceError{};
};

struct BSplineKnownValues
{
    double tolerance{};
    Vector<double> controlPoints;
    Vector<double> knots;
    Vector<double> x;
    Vector<double> y;
};

struct BlackKnownValue
{
    double tolerance{};
    double t{};
    double discountFactor{};
    double forward{};
    double volatility{};
    double strike{};
    double price{};
};

ExamplePipeline readExamplePipeline(const std::filesystem::path& path);
SyntheticSVICalibration readSyntheticSVICalibration(const std::filesystem::path& path);
BSplineKnownValues readBSplineKnownValues(const std::filesystem::path& path);
BlackKnownValue readBlackKnownValue(const std::filesystem::path& path);
} // namespace uv::tests::golden
