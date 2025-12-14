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

namespace uv::math::interp
{
	/**
	 * @brief Piecewise cubic Hermite interpolation (PCHIP).
	 *
	 * Computes a monotone piecewise cubic Hermite interpolant using function values
	 * `ys` at strictly increasing knot locations `xs`.
	 *
	 * Derivatives at the knots are computed internally using a monotonicity‐
	 * preserving formula. The result is then evaluated at the query points `x`
	 * using the generic Hermite spline routine.
	 *
	 * ### Behavior
	 *
	 * - For x < xs.front(), the value ys.front() is returned (flat extrapolation).
	 * - For x > xs.back(), the value ys.back() is returned (flat extrapolation).
	 * - Within the domain, each interval [xs[i], xs[i+1]] is a cubic polynomial.
	 *
	 * ### Requirements
	 *
	 * - `xs` must be strictly increasing.
	 * - `xs.size() == ys.size() >= 2`.
	 *
	 * ### Complexity
	 *
	 * - Slope computation: O(N) for N = xs.size().
	 * - Evaluation: O(M log N) for M = x.size() using binary search.
	 *
	 * ### References
	 *
	 * - Fritsch, C. and Carlson, R. (1980). "Monotone Piecewise Cubic Interpolation".
	 * - SciPy: `scipy.interpolate.PchipInterpolator`
	 *
	 * @tparam T     Floating‐point type satisfying `std::floating_point`.
	 *
	 * @param x      Query points at which to evaluate the interpolant.
	 * @param xs     Strictly increasing knot locations.
	 * @param ys     Function values corresponding to `xs`.
	 *
	 * @return       Interpolated values at each point in `x`.
	 */
	template <std::floating_point T>
	Vector<T> pchipInterp(
		const Vector<T>& x,
		const Vector<T>& xs,
		const Vector<T>& ys
	);

	///**
	// * @brief 2D tensor-product PCHIP interpolation on a rectangular grid.
	// *
	// * Performs shape-preserving piecewise cubic Hermite interpolation (PCHIP)
	// * along both dimensions of a 2D dataset defined on a rectangular grid.
	// *
	// * ### Method
	// *
	// * 1. Interpolate along x-direction for each fixed y (row pass).
	// * 2. Transpose intermediate matrix so columns become rows.
	// * 3. Interpolate along y-direction using identical row logic.
	// * 4. Transpose back to obtain final layout (y.size() × x.size()).
	// *	 
	// * @tparam T    Floating-point type satisfying std::floating_point.
	// *
	// * @param x     Target grid along x-axis (e.g. strikes), size Mx.
	// * @param y     Target grid along y-axis (e.g. maturities), size My.
	// * @param xs    Source x-grid, strictly increasing, size Nx.
	// * @param ys    Source y-grid, strictly increasing, size Ny.
	// * @param zs    Rectangular matrix of samples f(ys[i], xs[j]).
	// *
	// * @return      Matrix of size y.size() × x.size() containing
	// *              2D PCHIP-interpolated values.
	// */
	//template <std::floating_point T>
	//Matrix<T> pchipInterp2D(const Vector<T>& x,
	//	const Vector<T>& y,
	//	const Vector<T>& xs,
	//	const Vector<T>& ys,
	//	const Matrix<T>& zs);

	/**
	 * @brief Cubic Hermite interpolation with precomputed tangent slopes.
	 *
	 * Computes a piecewise cubic Hermite interpolant using values (xs, ys) and
	 * user-supplied first derivatives (dydx) at each node. The routine returns a
	 * C1-continuous interpolated value at arbitrary query points x by constructing,
	 * on each interval [xs[i], xs[i+1]], the standard cubic polynomial:
	 *
	 *     H(x) = y[i]
	 *            + d[i] * dx
	 *            + c2[i] * (dx)²
	 *            + c3[i] * (dx)³
	 *
	 * where dx = (x − xs[i]) and the coefficients c₂ and c₃ are computed locally
	 * from node spacing and slope information.
	 *
	 *
	 * ### Requirements
	 *
	 * - xs must be strictly increasing.
	 * - xs.size() == ys.size() == dydx.size().
	 * - xs.size() >= 1 (degenerate cases xs.size() == 0 or 1 are handled specially).
	 *
	 * ### Behavior
	 *
	 * - For x < xs.front(), flat extrapolation.
	 * - For x > xs.back(), flat extrapolation.
	 *
	 * ### Complexity
	 *
	 * - Preprocessing over the knots (xs, ys, dydx) is O(N).
	 * - Each query x[i] is evaluated in O(log N) time using binary search.
	 *   Overall cost is O(N + M log N) for M = x.size().
	 *
	 * ### References
	 *
	 *   - "Monotone cubic interpolation"
	 *     https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
	 *
	 * @tparam T   Floating-point type satisfying std::floating_point.
	 *
	 * @param x     Query points at which interpolation should be evaluated.
	 * @param xs    Strictly increasing interpolation x-grid.
	 * @param ys    Function values corresponding to xs.
	 * @param dydx  First derivative estimates at each xs.
	 *
	 * @return      Vector of interpolated values at each entry in x using
	 *              cubic Hermite interpolation.
	 */
	template <std::floating_point T>
	Vector<T> hermiteSplineInterp(
		const Vector<T>& x,
		const Vector<T>& xs,
		const Vector<T>& ys,
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

	namespace details
	{	
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

		/**
		 * @brief Row-wise PCHIP interpolation of a 2D matrix.
		 *
		 * Performs 1D piecewise cubic Hermite interpolation (PCHIP) independently on
		 * each row of a 2D matrix. For each fixed row i, the routine interpolates
		 * entries on the source grid (xs, ys[i]) to a new grid x, preserving the
		 * number of rows and replacing the columns.
		 *
		 * Given a matrix:
		 *
		 *     ys[i][j] = f( xs[j] ),   i = 0..(nRows-1),  j = 0..(nXS-1)
		 *
		 * the output matrix has the form:
		 *
		 *     result[i][k] = f_interp( x[k] ),  k = 0..(x.size()-1)
		 *
		 * @tparam T    Floating-point type satisfying std::floating_point.
		 *
		 * @param x     Target x-grid for interpolation (size M).
		 * @param xs    Source x-grid (size C), must be strictly increasing.
		 * @param ys    Matrix of values indexed as ys[row][col].
		 *
		 * @return      Matrix with the same number of rows as ys and
		 *              x.size() columns containing interpolated values.
		 */
		template <std::floating_point T>
		Matrix<T> pchipInterpRows(const Vector<T>& x,
			const Vector<T>& xs,
			const Matrix<T>& ys);

	} // namespace details

} // namespace uv::math::interp

#include "Interpolation.inl"