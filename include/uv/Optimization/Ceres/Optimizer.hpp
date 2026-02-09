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

#include "Base/Types.hpp"
#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Policy.hpp"

#include <ceres/ceres.h>
#include <optional>

namespace uv::opt::ceres
{

template <typename PolicyT = Policy<>> class Optimizer
{
  private:
    Config config_;
    std::optional<Vector<double>> lowerBounds_;
    std::optional<Vector<double>> upperBounds_;
    Vector<double> x_;
    ::ceres::Problem problem_;

    bool isInitialized_{false};
    bool isRunStarted_{false};

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