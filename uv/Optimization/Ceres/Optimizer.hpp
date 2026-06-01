// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Policy.hpp"

#include <ceres/ceres.h>
#include <memory>
#include <optional>

namespace uv::opt::ceres
{

template <typename PolicyT = Policy<>> class Optimizer
{
  private:
    Config config_;
    std::unique_ptr<::ceres::LossFunction> loss_;
    ::ceres::Solver::Options options_;

    bool isInitialized_{false};
    bool isRunStarted_{false};

    std::optional<Vector<double>> lowerBounds_;
    std::optional<Vector<double>> upperBounds_;
    Vector<double> x_;
    ::ceres::Problem problem_;

    void setOptions_();

    void setBounds_(
        std::size_t n,
        std::span<const double> lowerBounds,
        std::span<const double> upperBounds
    );

    void clampStoredBounds_();

    void requireInitialized_() const;
    void requireRunStarted_() const;

  public:
    Optimizer() = delete;

    explicit Optimizer(const Config& config);

    void initialize(
        std::span<const double> initGuess,
        std::span<const double> lowerBounds,
        std::span<const double> upperBounds
    );

    void beginRun();

    void addResidualBlock(std::unique_ptr<::ceres::CostFunction> cf);

    void solveInPlace();

    std::span<const double> solve();

    std::span<const double> params() const;
};
} // namespace uv::opt::ceres

#include "Optimization/Ceres/Detail/Optimizer.inl"