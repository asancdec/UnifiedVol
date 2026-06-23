
// SPDX-License-Identifier: Apache-2.0

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

    std::array<double, N> lowerBounds_{};
    std::array<double, N> upperBounds_{};
    std::array<double, N> initGuess_{};

    NloptFunction userFn_{nullptr};
    void* userData_{nullptr};
    unsigned iterCount_{0U};

    std::optional<double> userValue_;

    [[gnu::hot]] static double
    objectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept;

  public:
    Optimizer() = delete;

    explicit Optimizer(const Config<N>& config);

    Optimizer<N, Algo> fresh() const;

    void setGuessBounds(
        std::array<double, N> initGuess,
        std::array<double, N> lowerBounds,
        std::array<double, N> upperBounds
    );

    void addInequalityConstraint(NloptFunction c, void* data);

    void addInequalityMConstraint(std::size_t m, NloptMFunction c, void* data);

    void setMinObjective(NloptFunction f, void* data);

    Vector<double> optimize();

    void setUserValue(double v) noexcept;

    const double& eps() const noexcept;

    double tol() const noexcept;

    const double& userValue() const;
};
} // namespace uv::opt::nlopt

#include "Optimization/NLopt/Detail/Optimizer.inl"
