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
    using enum LinearSolver;

    switch (a)
    {
    case DenseQR:
        return ::ceres::DENSE_QR;
    case DenseNormalCholesky:
        return ::ceres::DENSE_NORMAL_CHOLESKY;
    case SparseNormalCholesky:
        return ::ceres::SPARSE_NORMAL_CHOLESKY;
    case SparseSchur:
        return ::ceres::SPARSE_SCHUR;
    }
    UNREACHABLE(LinearSolver, a);
}

std::unique_ptr<::ceres::LossFunction> makeLoss(Loss a, double lossParam)
{
    using enum Loss;

    switch (a)
    {
    case None:
        return nullptr;
    case Huber:
        return std::make_unique<::ceres::HuberLoss>(lossParam);
    case Cauchy:
        return std::make_unique<::ceres::CauchyLoss>(lossParam);
    case SoftL1:
        return std::make_unique<::ceres::SoftLOneLoss>(lossParam);
    }

    UNREACHABLE(Loss, a);
}

} // namespace uv::opt::ceres::detail
