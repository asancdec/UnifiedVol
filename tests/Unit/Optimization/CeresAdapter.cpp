// SPDX-License-Identifier: Apache-2.0

#include "Optimization/Ceres/Detail/CeresAdapter.hpp"
#include "Base/Errors/Errors.hpp"

#include <ceres/ceres.h>
#include <gtest/gtest.h>

namespace ceres_detail = uv::opt::ceres::detail;
namespace ceres_opt = uv::opt::ceres;

TEST(UnitOptimizationCeresAdapter, MapsTrustRegionStrategies)
{
    EXPECT_EQ(
        ceres_detail::toCeres(ceres_opt::TrustRegionStrategy::LevenbergMarquardt),
        ::ceres::LEVENBERG_MARQUARDT
    );
    EXPECT_EQ(
        ceres_detail::toCeres(ceres_opt::TrustRegionStrategy::Dogleg),
        ::ceres::DOGLEG
    );
}

TEST(UnitOptimizationCeresAdapter, MapsLinearSolvers)
{
    EXPECT_EQ(ceres_detail::toCeres(ceres_opt::LinearSolver::DenseQR), ::ceres::DENSE_QR);
    EXPECT_EQ(
        ceres_detail::toCeres(ceres_opt::LinearSolver::DenseNormalCholesky),
        ::ceres::DENSE_NORMAL_CHOLESKY
    );
    EXPECT_EQ(
        ceres_detail::toCeres(ceres_opt::LinearSolver::SparseNormalCholesky),
        ::ceres::SPARSE_NORMAL_CHOLESKY
    );
    EXPECT_EQ(
        ceres_detail::toCeres(ceres_opt::LinearSolver::SparseSchur),
        ::ceres::SPARSE_SCHUR
    );
}

TEST(UnitOptimizationCeresAdapter, BuildsLossFunctions)
{
    EXPECT_EQ(ceres_detail::makeLoss(ceres_opt::Loss::None, 1.0), nullptr);
    EXPECT_NE(ceres_detail::makeLoss(ceres_opt::Loss::Huber, 1.0), nullptr);
    EXPECT_NE(ceres_detail::makeLoss(ceres_opt::Loss::Cauchy, 1.0), nullptr);
    EXPECT_NE(ceres_detail::makeLoss(ceres_opt::Loss::SoftL1, 1.0), nullptr);
}

TEST(UnitOptimizationCeresAdapter, RejectsUnknownEnumValues)
{
    EXPECT_THROW(
        ceres_detail::toCeres(static_cast<ceres_opt::TrustRegionStrategy>(99)),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        ceres_detail::toCeres(static_cast<ceres_opt::LinearSolver>(99)),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(ceres_detail::makeLoss(static_cast<ceres_opt::Loss>(99), 1.0)),
        uv::errors::UnifiedVolError
    );
}
