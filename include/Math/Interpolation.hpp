// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Interpolation.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
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

#include "Utils/Types.hpp"

#include <concepts>

namespace uv::math
{
	/**
	 * @brief Cubic Hermite interpolation with precomputed tangent slopes.
	 *
	 * Computes a piecewise cubic Hermite interpolant using values (xs, ys) and
	 * user-supplied first derivatives (dydx) at each node. The routine returns a
	 * C1-continuous interpolated value at arbitrary x by constructing, on each
	 * interval [x[i], x[i+1]], the standard cubic polynomial:
	 *
	 *     H(x) = y[i]
	 *            + d[i] * Δx
	 *            + c2[i] * (Δx)²
	 *            + c3[i] * (Δx)³
	 *
	 * where Δx = (x − x[i]) and the coefficients c₂ and c₃ are computed locally
	 * from node spacing and slope information.
	 *
	 * ### Features
	 *
	 * - Local: Each interval is computed independently, without global solves.
	 * - C¹ continuity: Function value and first derivative are continuous across
	 *   all spline knots.
	 * - Handles arbitrary (not necessarily monotone) derivative data.
	 * - Supports clamped extrapolation outside the domain [xs.front(), xs.back()].
	 *
	 * ### Requirements
	 *
	 * - xs must be strictly increasing.
	 * - xs.size() == ys.size() == dydx.size().
	 * - xs.size() >= 2.
	 *
	 * ### Behavior
	 *
	 * - For x < xs.front(), the result is the constant value ys.front().
	 * - For x > xs.back(), the result is the constant value ys.back().
	 * - No new memory allocations besides O(N) temporary buffers for interval data.
	 *
	 * ### Complexity
	 *
	 * Time complexity is O(N) for preprocessing (step sizes, secant slopes, and
	 * cubic coefficients) plus O(log N) for interval lookup using binary search.
	 *
	 * ### References
	 *
	 *   - "Monotone cubic interpolation"
	 *     https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
	 *
	 * @tparam T   Floating-point type satisfying std::floating_point.
	 *
	 * @param x     Query point at which interpolation should be evaluated.
	 * @param xs    Strictly increasing interpolation x-grid.
	 * @param ys    Function values corresponding to xs.
	 * @param dydx  First derivative estimates at each xs.
	 *
	 * @return      Interpolated value at x using cubic Hermite interpolation.
	 */
	template <std::floating_point T>
	T interpolateCubicHermiteSpline(const T xs,
		const Vector<T>& x,
		const Vector<T>& y,
		const Vector<T>& dydx
		);

	/**
	 * @brief Compute monotone PCHIP derivative estimates (Fritsch–Butland method).
	 *
	 * Produces first-derivative values suitable for monotone piecewise cubic Hermite
	 * interpolation (PCHIP). Given strictly increasing nodes x[i] and corresponding
	 * function values y[i], this routine computes a vector d[i] of tangent slopes
	 * such that a cubic Hermite interpolant constructed from (x[i], y[i], d[i]) is:
	 *
	 *   - Shape-preserving: no new extrema are introduced.
	 *   - Local: each d[i] depends only on neighboring data.
	 *   - One-pass and O(N).
	 *
	 * ### Mathematical overview
	 *
	 * Interior derivatives are computed using the weighted harmonic mean formula of
	 * Fritsch and Butland (1984), while endpoint derivatives use the standard
	 * three-point, one-sided shape-preserving scheme of Fritsch & Carlson (1980)
	 * as implemented in MATLAB's `pchip`.
	 *
	 * If secant slopes S[i] and S[i+1] change sign, the derivative d[i] is set to 0
	 * to prevent overshoot and enforce monotonicity.
	 *
	 * ### Special cases
	 *
	 * - If x.size() == 2, the result is simply a straight line segment and the
	 *   returned derivative vector is {S, S}, where S is the secant slope.
	 *
	 * ### Requirements
	 *
	 * - x must be strictly increasing.
	 * - x.size() == y.size().
	 * - x.size() >= 2.
	 *
	 * ### References
	 *
	 *   - F. N. Fritsch and J. Butland,
	 *     "A Method for Constructing Local Monotone Piecewise Cubic Interpolants,"
	 *     SIAM Journal on Scientific and Statistical Computing, 5(2), 300–304 (1984).
	 *
	 * @tparam T  Floating-point type.
	 * 
	 * @param xs   Strictly increasing interpolation grid.
	 * @param ys   Function values corresponding to x.
	 *
	 * @return    Vector d of tangent derivatives with size equal to x.size().
	 *
	 * @note Complexity: O(N), no global solves or matrix factorizations.
	 */
	template <std::floating_point T>
	Vector<T> pchipDerivatives(const Vector<T>& xs,
		const Vector<T>& ys);

	/**
	 * @brief Compute a shape-preserving PCHIP endpoint slope (one-sided derivative).
	 *
	 * This implements the non-centered, shape-preserving three-point endpoint formula
	 * used in the Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).
	 *
	 * Given three consecutive points x(i), x(i+1), x(i+2) with spacings
	 *   h1 = x(i+1) - x(i),
	 *   h2 = x(i+2) - x(i+1),
	 * and secant slopes
	 *   S1 = (y(i+1) - y(i))     / h1,
	 *   S2 = (y(i+2) - y(i+1))   / h2,
	 * this function returns a one-sided derivative d at the endpoint that:
	 *
	 *  - Uses the higher-order formula
	 *        d = ((2*h1 + h2)*S1 - h1*S2) / (h1 + h2),
	 *    which is O(h^2) accurate on a uniform grid.
	 *  - Enforces shape preservation / monotonicity by:
	 *        * Setting d = 0 if d has opposite sign to S1.
	 *        * Clamping |d| to 3*|S1| if S1 and S2 have opposite signs and
	 *          |d| > 3*|S1|.
	 *
	 * ### References
	 *
	 *   - F.N. Fritsch and R.E. Carlson,
	 *     "Monotone Piecewise Cubic Interpolation",
	 *     SIAM Journal on Numerical Analysis, 17(2), 1980.
	 *
	 *   - C. Moler, "Numerical Computing with MATLAB", SIAM, 2004.
	 *     DOI:10.1137/1.9780898717952
	 *
	 * @tparam T  Floating-point type.
	 * 
	 * @param h1  First interval length (x(i+1) - x(i)), must be > 0.
	 * @param h2  Second interval length (x(i+2) - x(i+1)), must be > 0.
	 * @param S1  First secant slope over [x(i),   x(i+1)].
	 * @param S2  Second secant slope over [x(i+1), x(i+2)].
	 * 
	 * @return    Shape-preserving one-sided derivative at the endpoint.
	 */
	template <std::floating_point T>
	T pchipEndpointSlope(const T h1,
		const T h2,
		const T S1,
		const T S2
	) noexcept;

} // namespace uv

#include "Interpolation.inl"