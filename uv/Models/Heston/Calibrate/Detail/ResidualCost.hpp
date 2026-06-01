// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/Heston/Calibrate/Detail/MaturitySlice.hpp"
#include "Models/Heston/Price/Pricer.hpp"
#include "Optimization/Ceres/Config.hpp"

#include <concepts>
#include <cstddef>
#include <memory>

#include <ceres/cost_function.h>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N> struct SliceCommon
{
    const MaturitySlice* s_{nullptr};
    const price::Pricer<T, N>* pricer_{nullptr};

    SliceCommon(const MaturitySlice& s, const price::Pricer<T, N>& pricer) noexcept;

    void residualOnly(const double* p, double* residuals) const;

    void residualAndJac(const double* p, double* residuals, double* J) const;
};

template <std::floating_point T, std::size_t N> class SliceJacobian final
    : public ::ceres::CostFunction
{
  public:
    explicit SliceJacobian(SliceCommon<T, N> c) noexcept;

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)
        const override;

  private:
    SliceCommon<T, N> common_;
};

template <std::floating_point T, std::size_t N> struct SliceFunctor
{
    SliceCommon<T, N> common_;

    bool operator()(double const* const* parameters, double* residuals) const;
};

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostAnalytic(const MaturitySlice& slice, const price::Pricer<T, N>& pricer);

template <std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeSliceCostNumericForward(
    const MaturitySlice& slice,
    const price::Pricer<T, N>& pricer
);

template <std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeSliceCostNumericCentral(
    const MaturitySlice& slice,
    const price::Pricer<T, N>& pricer
);

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction>
makeSliceCost(const MaturitySlice& slice, const price::Pricer<T, N>& pricer);
} // namespace uv::models::heston::calibrate::detail

#include "Models/Heston/Calibrate/Detail/ResidualCost.inl"