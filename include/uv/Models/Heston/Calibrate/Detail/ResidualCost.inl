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

#include "Base/Macros/Unreachable.hpp"

#include <ceres/dynamic_numeric_diff_cost_function.h>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N> SliceCommon<T, N>::SliceCommon(
    const MaturitySlice& s,
    const price::Pricer<T, N>& pricer
) noexcept
    : s_(&s),
      pricer_(&pricer)
{
}

template <std::floating_point T, std::size_t N>
void SliceCommon<T, N>::residualOnly(const double* p, double* residuals) const
{
    const std::size_t M = s_->K.size();

    for (std::size_t i = 0; i < M; ++i)
    {
        const double model =
            pricer_
                ->callPrice(p[0], p[1], p[2], p[3], p[4], s_->t, s_->dF, s_->F, s_->K[i]);

        residuals[i] = (model - s_->mkt[i]) * s_->w[i];
    }
}

template <std::floating_point T, std::size_t N> void
SliceCommon<T, N>::residualAndJac(const double* p, double* residuals, double* J) const
{
    const std::size_t M = s_->K.size();

    for (std::size_t i = 0; i < M; ++i)
    {
        const auto pg = pricer_->callPriceWithGradient(
            p[0],
            p[1],
            p[2],
            p[3],
            p[4],
            s_->t,
            s_->dF,
            s_->F,
            s_->K[i]
        );

        const double wi = s_->w[i];

        residuals[i] = (pg[0] - s_->mkt[i]) * wi;

        double* row = &J[i * 5];
        row[0] = pg[1] * wi;
        row[1] = pg[2] * wi;
        row[2] = pg[3] * wi;
        row[3] = pg[4] * wi;
        row[4] = pg[5] * wi;
    }
}

template <std::floating_point T, std::size_t N>
SliceJacobian<T, N>::SliceJacobian(SliceCommon<T, N> c) noexcept
    : common_(std::move(c))
{
    set_num_residuals(static_cast<int>(common_.s_->K.size()));
    mutable_parameter_block_sizes()->push_back(5);
}

template <std::floating_point T, std::size_t N> bool SliceJacobian<T, N>::Evaluate(
    double const* const* parameters,
    double* residuals,
    double** jacobians
) const
{
    const double* p = parameters[0];

    if (jacobians && jacobians[0])
    {
        common_.residualAndJac(p, residuals, jacobians[0]);
    }
    else
    {
        common_.residualOnly(p, residuals);
    }

    return true;
}

template <std::floating_point T, std::size_t N> bool
SliceFunctor<T, N>::operator()(double const* const* parameters, double* residuals) const
{
    common_.residualOnly(parameters[0], residuals);
    return true;
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostAnalytic(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    return std::make_unique<SliceJacobian<T, N>>(SliceCommon<T, N>{slice, pricer});
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostNumericForward(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    using Functor = SliceFunctor<T, N>;
    using Cost = ::ceres::DynamicNumericDiffCostFunction<Functor, ::ceres::FORWARD>;

    auto* cost = new Cost(new Functor{SliceCommon<T, N>{slice, pricer}});
    cost->AddParameterBlock(5);
    cost->SetNumResiduals(static_cast<int>(slice.K.size()));
    return std::unique_ptr<::ceres::CostFunction>(cost);
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostNumericCentral(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    using Functor = SliceFunctor<T, N>;
    using Cost = ::ceres::DynamicNumericDiffCostFunction<Functor, ::ceres::CENTRAL>;

    auto* cost = new Cost(new Functor{SliceCommon<T, N>{slice, pricer}});
    cost->AddParameterBlock(5);
    cost->SetNumResiduals(static_cast<int>(slice.K.size()));
    return std::unique_ptr<::ceres::CostFunction>(cost);
}

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction>
makeSliceCost(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    if constexpr (Mode == opt::ceres::GradientMode::Analytic)
    {
        return makeSliceCostAnalytic<T, N>(slice, pricer);
    }
    else if constexpr (Mode == opt::ceres::GradientMode::NumericForward)
    {
        return makeSliceCostNumericForward<T, N>(slice, pricer);
    }
    else if constexpr (Mode == opt::ceres::GradientMode::NumericCentral)
    {
        return makeSliceCostNumericCentral<T, N>(slice, pricer);
    }
    else [[unlikely]]
    {
        UV_UNREACHABLE(opt::ceres::GradientMode, Mode);
    }
}

} // namespace uv::models::heston::calibrate::detail