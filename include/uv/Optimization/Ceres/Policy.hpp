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

#include <ceres/ceres.h>
#include <memory>
#include <type_traits>

namespace uv::opt::ceres
{

template <
    typename LossType = void,
    ::ceres::TrustRegionStrategyType TrustRegionStrategy = ::ceres::LEVENBERG_MARQUARDT,
    ::ceres::LinearSolverType LinearSolver = ::ceres::DENSE_QR>
struct Policy
{
    using loss_type = LossType;

    static constexpr ::ceres::TrustRegionStrategyType trustRegionStrategy =
        TrustRegionStrategy;
    static constexpr ::ceres::LinearSolverType linearSolver = LinearSolver;

    static std::unique_ptr<::ceres::LossFunction> makeLoss(double lossParam)
    {
        if constexpr (std::is_void_v<LossType>)
        {

            return nullptr;
        }
        else
        {

            return std::make_unique<LossType>(lossParam);
        }
    }
};
} // namespace uv::opt::ceres
