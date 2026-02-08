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
#include "Math/Interpolation/Policies.hpp"

#include <concepts>
#include <span>

namespace uv::math::interp
{

template <std::floating_point T, class DerivPolicy, class EvalPolicy> struct Interpolator
{

    using value_type = T;
    using derivatives_type = DerivPolicy;
    using evaluator_type = EvalPolicy;

    DerivPolicy deriv{};

    EvalPolicy eval{};

    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> y,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;
};

template <std::floating_point T>
using PchipInterpolator = Interpolator<T, PchipDerivatives<T>, HermiteEval<T>>;

} // namespace uv::math::interp

#include "Math/Interpolation/Detail/Interpolator.inl"