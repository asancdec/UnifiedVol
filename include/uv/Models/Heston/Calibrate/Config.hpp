// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "Base/Types.hpp"
#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Policy.hpp"
#include "Optimization/Cost.hpp"

#include <cstddef>
#include <string_view>

namespace uv::models::heston::calibrate
{

struct Config
{
    double tolerance{1e-10};
    unsigned int maxEval{10000};
    opt::ceres::Verbosity verbosity{opt::ceres::Verbosity::Summary};
    opt::cost::WeightATM<double> weightATM{.wATM = 8.0, .k0 = 0.3};
    int numThreads{-1};
};

inline constexpr opt::ceres::GradientMode HestonGradient =
    opt::ceres::GradientMode::Analytic;

inline constexpr std::size_t defaultNodes{200};

using HestonPolicy = opt::ceres::Policy<
    opt::ceres::TrustRegionStrategy::LevenbergMarquardt,
    opt::ceres::LinearSolver::DenseQR>;

namespace detail
{
inline Vector<std::string_view> paramNames{"kappa", "theta", "sigma", "rho", "v0"};
}

} // namespace uv::models::heston::calibrate