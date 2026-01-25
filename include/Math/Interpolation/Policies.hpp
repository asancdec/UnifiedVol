// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Policies.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Interpolation policies and definitions.
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

#include <span>
#include <concepts>

namespace uv::math::interp
{
    /**
     * @brief Concept: policy provides node-derivative computation.
     *
     * Requires a member function:
     * `void derivatives(xs, ys, dydx, doValidate)`.
     *
     * @tparam D Policy type.
     * @tparam T Value type.
     */
    template<class D, class T>
    concept HasDerivatives = requires
    (
        const D& d,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> dydx,
        bool doValidate
    )
    {
        { d.derivatives(xs, ys, dydx, doValidate) } -> std::same_as<void>;
    };

    /**
     * @brief Concept: policy provides spline evaluation given node derivatives.
     *
     * Requires a member function:
     * `void evaluate(x, xs, ys, dydx, y, doValidate)`.
     *
     * @tparam D Policy type.
     * @tparam T Value type.
     */
    template <class D, class T>
    concept HasEvaluate = requires
    (
        const D& d,
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate
    )
    {
        { d.evaluate(x, xs, ys, dydx, y, doValidate) } -> std::same_as<void>;
    };

    /**
     * @brief CRTP base for derivative-computation policies.
     *
     * Provides:
     * - `operator()(xs, ys)` returning a fresh derivative vector.
     * - `operator()(xs, ys, dydx)` filling a provided buffer.
     *
     * @tparam Derived Concrete policy type (CRTP).
     * @tparam T Floating-point type.
     */
    template<class Derived, std::floating_point T>
    struct Derivatives
    {
        /**
         * @brief Compute node derivatives and return them as a vector.
         *
         * @param xs Knot locations (strictly increasing).
         * @param ys Knot values.
         * @param doValidate If true, validates inputs and throws on error.
         * @return Vector of size `xs.size()` with derivatives at knots.
         */
        Vector<T> operator()
        (
            std::span<const T> xs,
            std::span<const T> ys,
            bool doValidate = true
        ) const
            requires HasDerivatives<Derived, T>;

        /**
         * @brief Compute node derivatives into a caller-provided buffer.
         *
         * @param xs Knot locations (strictly increasing).
         * @param ys Knot values.
         * @param dydx Output buffer (size `xs.size()`).
         * @param doValidate If true, validates inputs and throws on error.
         */
        void operator()
        (
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<T> dydx,
            bool doValidate = true
        ) const
            requires HasDerivatives<Derived, T>;
    };

    /**
     * @brief PCHIP derivative policy (shape-preserving).
     *
     * Computes knot derivatives suitable for monotone/shape-preserving
     * cubic Hermite interpolation.
     *
     * @tparam T Floating-point type.
     */
    template <std::floating_point T>
    struct PchipDerivatives : Derivatives<PchipDerivatives<T>, T>
    {
        /**
         * @brief Compute PCHIP knot derivatives.
         *
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Output buffer for derivatives (size `xs.size()`).
         * @param doValidate If true, validates inputs and throws on error.
         */
        void derivatives
        (
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<T> dydx,
            bool doValidate
        ) const;
    };

    /**
     * @brief CRTP base for evaluation policies.
     *
     * Provides:
     * - vector evaluation `operator()(x, xs, ys, dydx)`
     * - in-place evaluation `operator()(x, xs, ys, dydx, y)`
     * - scalar evaluation `operator()(x0, xs, ys, dydx)`
     *
     * @tparam Derived Concrete policy type (CRTP).
     * @tparam T Floating-point type.
     */
    template<class Derived, std::floating_point T>
    struct Evaluate
    {
        /**
         * @brief Evaluate interpolant at multiple points (allocates output).
         *
         * @param x Query points.
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Knot derivatives.
         * @param doValidate If true, validates inputs and throws on error.
         * @return Values at `x` (size `x.size()`).
         */
        Vector<T> operator()
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            bool doValidate = true
        ) const
            requires HasEvaluate<Derived, T>;

        /**
         * @brief Evaluate interpolant at multiple points into a buffer.
         *
         * @param x Query points.
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Knot derivatives.
         * @param y Output buffer (size `x.size()`).
         * @param doValidate If true, validates inputs and throws on error.
         */
        void operator()
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            std::span<T> y,
            bool doValidate = true
        ) const
            requires HasEvaluate<Derived, T>;

        /**
         * @brief Evaluate interpolant at a single point.
         *
         * @param x Query point.
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Knot derivatives.
         * @param doValidate If true, validates inputs and throws on error.
         * @return Interpolated value at `x`.
         */
        T operator()
        (
            T x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            bool doValidate = true
        ) const
            requires HasEvaluate<Derived, T>;
    };

    /**
     * @brief Cubic Hermite spline evaluation policy.
     *
     * Evaluates a cubic Hermite interpolant given knots (xs, ys)
     * and knot derivatives (dydx).
     *
     * @tparam T Floating-point type.
     */
    template<std::floating_point T>
    struct HermiteEval : Evaluate<HermiteEval<T>, T>
    {
        /**
         * @brief Evaluate cubic Hermite spline at points x.
         *
         * @param x Query points.
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Knot derivatives.
         * @param y Output buffer (size `x.size()`).
         * @param doValidate If true, validates inputs and throws on error.
         */
        void evaluate
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            std::span<T> y,
            bool doValidate
        ) const;
    };

    namespace detail
    {
        /**
         * @brief Core cubic Hermite spline implementation.
         *
         * Handles validation (optional), interval search, interpolation, and
         * out-of-bounds behavior (current implementation uses flat extrapolation).
         *
         * @tparam T Floating-point type.
         */
        template <std::floating_point T>
        void hermiteSplineInterp
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            std::span<T> y,
            bool doValidate
        );

        /**
         * @brief Compute shape-preserving PCHIP knot derivatives.
         *
         * @param xs Knot locations (strictly increasing).
         * @param ys Knot values.
         * @param dydx Output buffer (size `xs.size()`).
         * @param doValidate If true, validates inputs and throws on error.
         */
        template <std::floating_point T>
        void pchipDerivatives
        (
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<T> dydx,
            bool doValidate
        );

        /**
         * @brief Shape-preserving one-sided endpoint derivative for PCHIP.
         *
         * @param h1 First interval length (> 0).
         * @param h2 Second interval length (> 0).
         * @param S1 Secant slope over first interval.
         * @param S2 Secant slope over second interval.
         * @return Endpoint derivative.
         */
        template <std::floating_point T>
        T pchipEndpointSlope
        (
            T h1,
            T h2,
            T S1,
            T S2
        ) noexcept;

        /**
         * @brief Validate inputs for derivative computation.
         *
         * Checks sizes, finiteness, and strict monotonicity of xs.
         *
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Derivative buffer (size check only).
         */
        template <std::floating_point T>
        void validateInputsDerivatives
        (
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx
        );

        /**
         * @brief Validate inputs for spline evaluation.
         *
         * Validates knot grid and derivative array, checks x/y sizes, and
         * ensures values are finite.
         *
         * @param x Query points.
         * @param xs Knot locations.
         * @param ys Knot values.
         * @param dydx Knot derivatives.
         * @param y Output buffer (size check only).
         */
        template <std::floating_point T>
        void validateInputsEvaluate
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            std::span<const T> dydx,
            std::span<const T> y
        );

    } // namespace detail
} // namespace uv::math::interp

#include "Policies.inl"