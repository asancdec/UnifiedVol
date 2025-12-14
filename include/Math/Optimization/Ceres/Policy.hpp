// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Policy.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

namespace uv::math::opt::ceres
{   
    /**
     * @brief Compile-time policy for configuring Ceres optimisation behaviour.
     *
     * Encapsulates solver strategy selection and optional robust loss
     * construction via template parameters.
     *
     * @tparam LossType Robust loss type (e.g. ceres::HuberLoss). Use `void` for no loss.
     * @tparam TrustRegionStrategy Ceres trust-region strategy (default: Levenberg–Marquardt).
     * @tparam LinearSolver Linear solver type used by Ceres.
     */
    template    
    <
    typename LossType = void,
    ::ceres::TrustRegionStrategyType TrustRegionStrategy = ::ceres::LEVENBERG_MARQUARDT,
    ::ceres::LinearSolverType LinearSolver = ::ceres::DENSE_QR
    >
    struct Policy
    {
        // ---------- Solver configurations ----------

        static constexpr ::ceres::TrustRegionStrategyType trustRegionStrategy = TrustRegionStrategy;
        static constexpr ::ceres::LinearSolverType linearSolver = LinearSolver;

        // ---------- Robust loss construction ----------

        static std::unique_ptr<::ceres::LossFunction> makeLoss(double lossParam)
        {
            if constexpr (std::is_same_v<LossType, void>) 
            {
                // No robust loss
                return nullptr; 
            }
            else 
            {
                // Construct loss with scale parameter
                return std::make_unique<LossType>(lossParam);
            }
        }
    };
}

