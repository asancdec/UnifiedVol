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
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#pragma once

#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Detail/CeresAdapter.hpp"

#include <ceres/ceres.h>
#include <memory>

namespace uv::opt::ceres
{
template <
    TrustRegionStrategy TR = TrustRegionStrategy::LevenbergMarquardt,
    LinearSolver LS = LinearSolver::DenseQR,
    Loss L = Loss::None>
struct Policy
{
    inline static ::ceres::TrustRegionStrategyType trustRegionStrategy =
        ::uv::opt::ceres::detail::toCeres(TR);

    inline static ::ceres::LinearSolverType linearSolver =
        ::uv::opt::ceres::detail::toCeres(LS);

    static std::unique_ptr<::ceres::LossFunction> makeLoss(double lossParam)
    {
        return ::uv::opt::ceres::detail::makeLoss(L, lossParam);
    }
};
} // namespace uv::opt::ceres