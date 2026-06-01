// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Base/Types.hpp>

#include <string_view>

namespace uv::opt::ceres
{

enum class TrustRegionStrategy
{
    LevenbergMarquardt,
    Dogleg
};

enum class LinearSolver
{
    DenseQR,
    DenseNormalCholesky,
    SparseNormalCholesky,
    SparseSchur
};

enum class Loss
{
    None,
    Huber,
    Cauchy,
    SoftL1
};

enum class GradientMode
{
    Analytic,
    NumericForward,
    NumericCentral
};

enum class Verbosity
{
    None,
    Summary,
    FullReport
};

struct Config
{

    unsigned maxEval;
    double functionTol;
    double paramTol;
    double gradientTol;

    double lossScale;

    Vector<std::string_view> paramNames;

    Verbosity verbosity;

    int numThreads{1};
};

} // namespace uv::opt::ceres