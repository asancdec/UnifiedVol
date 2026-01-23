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
	 * @brief Policy for generating monotone PCHIP first derivatives.
	 *
	 * Guarantees that all users of this policy compute node derivatives
	 * using the same PCHIP scheme. Implementation is delegated to
	 * detail:: namespace.
	 */
	template <std::floating_point T>
    struct PchipDerivatives
    {
		Vector<T> operator()
		(
			std::span<const T> xs,
			std::span<const T> ys,
			bool doValidate = true
	    ) const;
    };

	/**
	 * @brief Policy for evaluating cubic Hermite interpolation.
	 *
	 * Guarantees consistent Hermite interpolation given precomputed
	 * node derivatives. Actual evaluation logic is implemented in
	 * the detail:: namespace.
	 */
    template <std::floating_point T>
    struct HermiteInterpolator
    {
		Vector<T> operator()
		(
			std::span<const T> x,
			std::span<const T> xs,
			std::span<const T> ys,
			std::span<const T> dydx,
			bool doValidate = true
		) const;
    };

	/**
	 * @brief Policy for PCHIP interpolation.
	 *
	 * Combines PCHIP derivative construction and Hermite evaluation to
	 * ensure a consistent, shape-preserving interpolation scheme across
	 * all consumers. The numerical details are implemented in detail:: 
	 * namespace.
	 */
    template <std::floating_point T>
    struct PchipInterpolator
    {
        Vector<T> operator()
        (
            std::span<const T> x,
            std::span<const T> xs,
            std::span<const T> ys,
            bool doValidate = true
        ) const;

        T operator()
        (
            T x,
            std::span<const T> xs,
            std::span<const T> ys,
            bool doValidate = true
        ) const;
    };

	namespace detail
	{	
		/**
		 * @brief Piecewise cubic Hermite interpolation (PCHIP) for multiple query points.
		 *
		 * Computes PCHIP node derivatives from (xs, ys) and evaluates the resulting
		 * Hermite spline at query points x. Uses flat extrapolation outside
		 * [xs.front(), xs.back()].
		 *
		 * Requirements:
		 * - xs is strictly increasing
		 * - xs.size() == ys.size() >= 2
		 *
		 * Complexity:
		 * - O(N) to compute derivatives, O(M log N) to evaluate (binary search)
		 *
		 * @tparam T Floating-point type.
		 * @param x Query points.
		 * @param xs Strictly increasing knot locations.
		 * @param ys Knot values corresponding to xs.
		 * @param doValidate If true, validate inputs and throw on invalid data.
		 * @return Interpolated values at each point in x.
		 */
		template <std::floating_point T>
		Vector<T> pchipInterp(std::span<const T> x,
			std::span<const T> xs,
			std::span<const T> ys,
			bool doValidate);

		/**
		 * @brief Piecewise cubic Hermite interpolation (PCHIP) for a single query point.
		 *
		 * Convenience overload of pchipInterp that evaluates at one x. Uses flat
		 * extrapolation outside [xs.front(), xs.back()].
		 *
		 * Requirements:
		 * - xs is strictly increasing
		 * - xs.size() == ys.size() >= 2
		 *
		 * Complexity:
		 * - O(N) to compute derivatives, O(log N) to evaluate
		 *
		 * @tparam T Floating-point type.
		 * @param x Query point.
		 * @param xs Strictly increasing knot locations.
		 * @param ys Knot values corresponding to xs.
		 * @param doValidate If true, validate inputs and throw on invalid data.
		 * @return Interpolated value at x.
		 */
		template <std::floating_point T>
		T pchipInterp(T x,
			std::span<const T> xs,
			std::span<const T> ys,
			bool doValidate);

		/**
		 * @brief Evaluate a cubic Hermite spline with precomputed node derivatives.
		 *
		 * Evaluates the piecewise cubic defined by knots (xs, ys) and node slopes dydx.
		 * Uses flat extrapolation outside [xs.front(), xs.back()].
		 *
		 * Requirements:
		 * - xs is strictly increasing
		 * - xs.size() == ys.size() == dydx.size() >= 2
		 * - all inputs are finite if validateInputs is true
		 *
		 * Complexity:
		 * - O(N) preprocessing, O(M log N) evaluation
		 *
		 * @tparam T Floating-point type.
		 * @param x Query points.
		 * @param xs Strictly increasing knot locations.
		 * @param ys Knot values corresponding to xs.
		 * @param dydx Node derivatives at xs.
		 * @param doValidate If true, validate inputs and throw on invalid data.
		 * @return Interpolated values at each point in x.
		 */
		template <std::floating_point T>
		Vector<T> hermiteSplineInterp(std::span<const T> x,
			std::span<const T> xs,
			std::span<const T> ys,
			std::span<const T> dydx,
			bool doValidate);

		/**
		 * @brief Compute shape-preserving PCHIP node derivatives.
		 *
		 * Computes first-derivative values at the knots (xs, ys) suitable for
		 * piecewise cubic Hermite interpolation (PCHIP). Interior derivatives
		 * are computed using a weighted harmonic mean when adjacent secant
		 * slopes have the same sign; otherwise the derivative is set to zero.
		 * Endpoint derivatives use a one-sided, shape-preserving rule.
		 *
		 * Special case:
		 * - If xs.size() == 2, returns {S, S}, where S is the secant slope.
		 *
		 * Requirements:
		 * - xs is strictly increasing
		 * - xs.size() == ys.size() >= 2
		 *
		 * @tparam T Floating-point type.
		 * @param xs Knot locations.
		 * @param ys Knot values.
		 * @param doValidate If true, validate inputs and throw on error.
		 * @return Vector of node derivatives with size xs.size().
		 *
		 * @note Complexity: O(N).
		 *
		 * References:
		 * - Fritsch and Carlson (1980), Monotone Piecewise Cubic Interpolation.
		 * - MATLAB PCHIP implementation.
		 */
		template <std::floating_point T>
		Vector<T> pchipDerivatives(std::span<const T> xs,
			std::span<const T> ys,
			bool doValidate);

		/**
		 * @brief Compute a shape-preserving PCHIP endpoint derivative.
		 *
		 * Computes a one-sided derivative at an endpoint using three consecutive
		 * points. The initial estimate is clamped to enforce monotonicity:
		 * - Set to 0 if it disagrees in sign with the adjacent secant.
		 * - Clamped to 3*|S1| if neighboring secants have opposite signs.
		 *
		 * @tparam T Floating-point type.
		 * @param h1 First interval length (must be > 0).
		 * @param h2 Second interval length (must be > 0).
		 * @param S1 Secant slope on the first interval.
		 * @param S2 Secant slope on the second interval.
		 * @return Shape-preserving endpoint derivative.
		 *
		 * References:
		 * - Fritsch and Carlson (1980).
		 * - Moler, Numerical Computing with MATLAB (2004).
		 */
		template <std::floating_point T>
		T pchipEndpointSlope(T h1,
			T h2,
			T S1,
			T S2) noexcept;

		/**
		 * @brief Validate interpolation grid inputs.
		 *
		 * Checks that xs and ys have equal size, contain at least two points,
		 * consist only of finite values, and that xs is strictly increasing.
		 *
		 * @tparam T Floating-point type.
		 * @param xs Interpolation grid.
		 * @param ys Values at grid points.
		 *
		 * @throws ErrorCode::InvalidArgument on invalid input.
		 */
		template <std::floating_point T>
		void validateInputs(std::span<const T> xs,
			std::span<const T> ys);

		/**
		 * @brief Validate inputs for Hermite spline interpolation.
		 *
		 * Extends validateInputs(xs, ys) by also checking that:
		 * - dydx has the same size as xs
		 * - x and dydx contain only finite values
		 *
		 * @tparam T Floating-point type.
		 * @param x Query points.
		 * @param xs Knot locations.
		 * @param ys Knot values.
		 * @param dydx Node derivatives.
		 *
		 * @throws ErrorCode::InvalidArgument on invalid input.
		 */
		template <std::floating_point T>
		void validateInputsHermiteSpline(std::span<const T> x,
			std::span<const T> xs,
			std::span<const T> ys,
			std::span<const T> dydx);

	} // namespace detail
} // namespace uv::math::interp

#include "Policies.inl"