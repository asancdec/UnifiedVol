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
#include "Optimization/Ceres/Config.hpp"

#include <ceres/ceres.h>
#include <concepts>
#include <memory>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N> struct ResidualCommon
{
    double t_;
    double dF_;
    double F_;
    double K_;
    double mkt_;
    double w_;
    const price::Pricer<T, N>* pricer_;

    ResidualCommon(
        double t,
        double dF,
        double F,
        double K,
        double callPriceMkt,
        double w,
        const price::Pricer<T, N>& pricer
    ) noexcept;

    double residualOnly(const double* p) const;

    void residualAndJac(const double* p, double* r, double* J) const;
};

template <std::floating_point T, std::size_t N> struct ResidualJacobian;

template <std::floating_point T, std::size_t N> struct ResidualFunctor;

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeCost(ResidualCommon<T, N> c);

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction> makeCost(
    double t,
    double dF,
    double F,
    double K,
    double callPriceMkt,
    double w,
    const price::Pricer<T, N>& pricer
);
} // namespace uv::models::heston::calibrate::detail

#include "Models/Heston/Calibrate/Detail/ResidualCost.inl"