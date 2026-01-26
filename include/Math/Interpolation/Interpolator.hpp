// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Interpolator.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Interpolator class header.
 *
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

#include "Core/Types.hpp"
#include "Math/Interpolation/Policies.hpp"

#include <concepts>
#include <span>

namespace uv::math::interp
{
/**
 * @brief Generic 1D interpolator with separable derivative and evaluation
 * policies.
 *
 * Combines:
 *  - a derivative policy (computes node derivatives)
 *  - an evaluation policy (evaluates the interpolant)
 *
 * Provides vector and scalar interpolation, with optional reuse
 * of precomputed derivatives.
 *
 * @tparam T           Floating-point type.
 * @tparam DerivPolicy Derivative computation policy.
 * @tparam EvalPolicy  Interpolation evaluation policy.
 */
template <std::floating_point T, class DerivPolicy, class EvalPolicy> struct Interpolator
{
    // ---------- Expose types ----------

    using value_type = T;
    using derivatives_type = DerivPolicy;
    using evaluator_type = EvalPolicy;

    // ---------- Members ----------

    /// Derivative computation policy
    DerivPolicy deriv{};

    /// Evaluation policy
    EvalPolicy eval{};

    // ---------- Core ----------

    /**
     * @brief Evaluate interpolant using precomputed derivatives.
     *
     * @param x           Query points.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param dydx        Node derivatives.
     * @param y           Output buffer (size = x.size()).
     * @param doValidate  Enable input validation.
     */
    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    // ---------- Convenience ----------

    /**
     * @brief Evaluate interpolant (derivatives computed internally).
     *
     * @param x           Query points.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param doValidate  Enable input validation.
     * @return Interpolated values at x.
     */
    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    /**
     * @brief Evaluate interpolant (derivatives computed internally) into an
     * output buffer.
     *
     * Computes PCHIP/Hermite (or whatever DerivPolicy/EvalPolicy are) node
     * derivatives internally, then evaluates the interpolant at @p x, writing
     * results into @p y.
     *
     * @param x           Query points.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param y           Output buffer (size = x.size()).
     * @param doValidate  Enable input validation.
     */
    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> y,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    /**
     * @brief Evaluate interpolant with externally provided derivatives.
     *
     * @param x           Query points.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param dydx        Node derivatives.
     * @param doValidate  Enable input validation.
     * @return Interpolated values at x.
     */
    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    /**
     * @brief Evaluate interpolant at a single point.
     *
     * Derivatives are computed internally.
     *
     * @param x           Query point.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param doValidate  Enable input validation.
     * @return Interpolated value.
     */
    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    /**
     * @brief Evaluate interpolant at a single point using given derivatives.
     *
     * @param x           Query point.
     * @param xs          Node locations.
     * @param ys          Node values.
     * @param dydx        Node derivatives.
     * @param doValidate  Enable input validation.
     * @return Interpolated value.
     */
    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;
};

/**
 * @brief PCHIP interpolator.
 *
 * Combines:
 *  - shape-preserving PCHIP derivatives
 *  - cubic Hermite spline evaluation
 */
template <std::floating_point T>
using PchipInterpolator = Interpolator<T, PchipDerivatives<T>, HermiteEval<T>>;

} // namespace uv::math::interp

#include "Interpolator.inl"