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

#include "Optimization/Ceres/Detail/CeresAdapter.hpp"

#include "Base/Macros/Unreachable.hpp"

namespace uv::opt::ceres::detail
{

::ceres::TrustRegionStrategyType toCeres(TrustRegionStrategy a)
{
    switch (a)
    {
    case TrustRegionStrategy::LevenbergMarquardt:
        return ::ceres::LEVENBERG_MARQUARDT;
    case TrustRegionStrategy::Dogleg:
        return ::ceres::DOGLEG;
    }

    UV_UNREACHABLE(TrustRegionStrategy, a);
}

::ceres::LinearSolverType toCeres(LinearSolver a)
{
    switch (a)
    {
    case LinearSolver::DenseQR:
        return ::ceres::DENSE_QR;
    case LinearSolver::DenseNormalCholesky:
        return ::ceres::DENSE_NORMAL_CHOLESKY;
    case LinearSolver::SparseNormalCholesky:
        return ::ceres::SPARSE_NORMAL_CHOLESKY;
    case LinearSolver::SparseSchur:
        return ::ceres::SPARSE_SCHUR;
    }
    UV_UNREACHABLE(LinearSolver, a);
}

std::unique_ptr<::ceres::LossFunction> makeLoss(Loss a, double lossParam)
{
    switch (a)
    {
    case Loss::None:
        return nullptr;
    case Loss::Huber:
        return std::make_unique<::ceres::HuberLoss>(lossParam);
    case Loss::Cauchy:
        return std::make_unique<::ceres::CauchyLoss>(lossParam);
    case Loss::SoftL1:
        return std::make_unique<::ceres::SoftLOneLoss>(lossParam);
    }

    UV_UNREACHABLE(Loss, a);
}

} // namespace uv::opt::ceres::detail