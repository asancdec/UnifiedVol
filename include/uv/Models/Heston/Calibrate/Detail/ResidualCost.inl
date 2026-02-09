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

#include <utility>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N>
ResidualCommon<T, N>::ResidualCommon(
    double t,
    double dF,
    double F,
    double K,
    double callPriceMkt,
    double w,
    const price::Pricer<T, N>& pricer
) noexcept
    : t_(t),
      dF_(dF),
      F_(F),
      K_(K),
      mkt_(callPriceMkt),
      w_(w),
      pricer_(&pricer)
{
}

template <std::floating_point T, std::size_t N>
double ResidualCommon<T, N>::residualOnly(const double* p) const
{
    const double model =
        pricer_->callPrice(p[0], p[1], p[2], p[3], p[4], t_, dF_, F_, K_);
    return (model - mkt_) * w_;
}

template <std::floating_point T, std::size_t N>
void ResidualCommon<T, N>::residualAndJac(const double* p, double* r, double* J) const
{
    const auto pg =
        pricer_->callPriceWithGradient(p[0], p[1], p[2], p[3], p[4], t_, dF_, F_, K_);

    r[0] = (pg[0] - mkt_) * w_;
    J[0] = pg[1] * w_;
    J[1] = pg[2] * w_;
    J[2] = pg[3] * w_;
    J[3] = pg[4] * w_;
    J[4] = pg[5] * w_;
}

template <std::floating_point T, std::size_t N>
struct ResidualJacobian final : public ::ceres::SizedCostFunction<1, 5>
{
    ResidualCommon<T, N> c_;

    explicit ResidualJacobian(ResidualCommon<T, N> c) noexcept
        : c_(std::move(c))
    {
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)
        const override
    {
        const double* p = parameters[0];

        if (jacobians && jacobians[0])
        {
            c_.residualAndJac(p, residuals, jacobians[0]);
        }
        else
        {
            residuals[0] = c_.residualOnly(p);
        }
        return true;
    }
};

template <std::floating_point T, std::size_t N> struct ResidualFunctor
{
    ResidualCommon<T, N> c_;

    bool operator()(double const* const p, double* residual) const
    {
        residual[0] = c_.residualOnly(p);
        return true;
    }
};

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeCost(ResidualCommon<T, N> c)
{
    if constexpr (Mode == opt::ceres::GradientMode::Analytic)
    {
        return std::make_unique<ResidualJacobian<T, N>>(std::move(c));
    }
    else if constexpr (Mode == opt::ceres::GradientMode::NumericForward)
    {
        using Functor = ResidualFunctor<T, N>;
        using Cost = ::ceres::NumericDiffCostFunction<Functor, ::ceres::FORWARD, 1, 5>;
        return std::unique_ptr<::ceres::CostFunction>(new Cost(new Functor{std::move(c)})
        );
    }
    else if constexpr (Mode == opt::ceres::GradientMode::NumericCentral)
    {
        using Functor = ResidualFunctor<T, N>;
        using Cost = ::ceres::NumericDiffCostFunction<Functor, ::ceres::CENTRAL, 1, 5>;
        return std::unique_ptr<::ceres::CostFunction>(new Cost(new Functor{std::move(c)})
        );
    }
    else [[unlikely]]
    {
        UV_UNREACHABLE(opt::ceres::GradientMode, Mode);
    }
}

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeCost(
    double t,
    double dF,
    double F,
    double K,
    double callPriceMkt,
    double w,
    const price::Pricer<T, N>& pricer
)
{
    return makeCost<Mode, T, N>(ResidualCommon<T, N>{t, dF, F, K, callPriceMkt, w, pricer}
    );
}

} // namespace uv::models::heston::calibrate::detail