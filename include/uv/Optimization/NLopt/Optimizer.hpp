
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
#include "Base/Utils/StopWatch.hpp"
#include "Optimization/NLopt/Algorithm.hpp"
#include "Optimization/NLopt/Config.hpp"

#include <array>
#include <cstddef>
#include <nlopt.hpp>
#include <optional>

namespace uv::opt::nlopt
{
template <std::size_t N, Algorithm Algo> class Optimizer
{
  private:
    using NloptFunction = double (*)(unsigned, const double*, double*, void*);
    using NloptMFunction = void (*)(
        unsigned m,
        double* result,
        unsigned n,
        const double* x,
        double* grad,
        void* data
    );

    Config<N> config_;
    ::nlopt::opt opt_;
    utils::StopWatch timer_;

    std::array<double, N> lowerBounds_;
    std::array<double, N> upperBounds_;
    std::array<double, N> initGuess_;

    NloptFunction userFn_;
    void* userData_;
    unsigned iterCount_;

    std::optional<double> userValue_;

    [[gnu::hot]] static double
    objectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept;

  public:
    Optimizer() = delete;

    explicit Optimizer(const Config<N>& config);

    Optimizer<N, Algo> fresh() const noexcept;

    void setGuessBounds(
        std::array<double, N> initGuess,
        std::array<double, N> lowerBounds,
        std::array<double, N> upperBounds
    ) noexcept;

    void addInequalityConstraint(NloptFunction c, void* data) noexcept;

    void addInequalityMConstraint(std::size_t m, NloptMFunction c, void* data) noexcept;

    void setMinObjective(NloptFunction f, void* data) noexcept;

    Vector<double> optimize();

    void setUserValue(double v) noexcept;

    const double& eps() const noexcept;

    double tol() const noexcept;

    const double& userValue() const noexcept;
};
} // namespace uv::opt::nlopt

#include "Optimization/NLopt/Detail/Optimizer.inl"