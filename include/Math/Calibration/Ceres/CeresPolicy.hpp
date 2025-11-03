/**
* CeresPolicy.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <ceres/ceres.h>
#include <memory>
#include <type_traits>

namespace uv
{
    template    
    <
    ::ceres::NumericDiffMethodType Method = ceres::CENTRAL,
    int M = 1,
    typename LossType = void,
    ::ceres::TrustRegionStrategyType TrustRegionStrategy = ::ceres::LEVENBERG_MARQUARDT,
    ::ceres::LinearSolverType LinearSolver = ::ceres::DENSE_QR
    >
    struct CeresPolicy
    {
        // Differentiation configuration
        static constexpr ::ceres::NumericDiffMethodType method = Method;
        static constexpr int m = M;

        // Solver configuration
        static constexpr ::ceres::TrustRegionStrategyType trustRegionStrategy = TrustRegionStrategy;
        static constexpr ::ceres::LinearSolverType linearSolver = LinearSolver;

        // Build the loss or return nullptr when Loss=void
        static std::unique_ptr<::ceres::LossFunction> makeLoss(double lossParam)
        {
            if constexpr (std::is_same_v<LossType, void>) 
            {
                return nullptr;
            }
            else 
            {
                return std::make_unique<LossType>(lossParam);
            }
        }
    };
}

