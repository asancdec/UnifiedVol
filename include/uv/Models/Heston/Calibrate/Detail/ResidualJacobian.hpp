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

#include "Models/Heston/Price/Pricer.hpp"
#include <ceres/ceres.h>
#include <concepts>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N>
struct ResidualJacobian final : public ceres::SizedCostFunction<1, 5>
{
    const double t_, dF_, F_, K_, callPriceMkt_, w_;
    const price::Pricer<T, N>* pricer_;

    ResidualJacobian(
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
          w_(w),
          callPriceMkt_(callPriceMkt),
          pricer_(&pricer)
    {
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)
        const override
    {
        const double* p = parameters[0];
        const auto pg =
            pricer_->callPriceWithGradient(p[0], p[1], p[2], p[3], p[4], t_, dF_, F_, K_);

        residuals[0] = double((pg[0] - callPriceMkt_) * w_);

        if ((jacobians != nullptr) && (jacobians[0] != nullptr))
        {
            double* J = jacobians[0];
            J[0] = double(pg[1] * w_);
            J[1] = double(pg[2] * w_);
            J[2] = double(pg[3] * w_);
            J[3] = double(pg[4] * w_);
            J[4] = double(pg[5] * w_);
        }
        return true;
    }
};
} // namespace uv::models::heston::calibrate::detail