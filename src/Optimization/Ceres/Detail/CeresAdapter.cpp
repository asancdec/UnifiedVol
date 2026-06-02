// SPDX-License-Identifier: Apache-2.0

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

    UNREACHABLE(TrustRegionStrategy, a);
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
    UNREACHABLE(LinearSolver, a);
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

    UNREACHABLE(Loss, a);
}

} // namespace uv::opt::ceres::detail