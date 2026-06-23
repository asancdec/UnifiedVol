// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Policy.hpp"
#include "Optimization/Cost.hpp"

#include <cstddef>

namespace uv::models::heston::calibrate
{

struct Config
{
    double tolerance{1e-11};
    unsigned int maxEval{10000};
    opt::ceres::Verbosity verbosity{opt::ceres::Verbosity::Summary};
    opt::cost::WeightATM<double> weightATM{};
    int numThreads{-1};
};

inline constexpr opt::ceres::GradientMode HestonGradient =
    opt::ceres::GradientMode::Analytic;

inline constexpr std::size_t defaultNodes{300};

using HestonPolicy = opt::ceres::Policy<
    opt::ceres::TrustRegionStrategy::LevenbergMarquardt,
    opt::ceres::LinearSolver::DenseQR>;

} // namespace uv::models::heston::calibrate
