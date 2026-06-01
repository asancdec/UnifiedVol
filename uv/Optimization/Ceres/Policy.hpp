// SPDX-License-Identifier: Apache-2.0

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