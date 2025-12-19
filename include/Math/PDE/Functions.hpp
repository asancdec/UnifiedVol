// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include "Utils/Aux/Errors.hpp"
#include "Math/PDE/TriDiag.hpp"
#include "Core/Matrix/Matrix.hpp"

#include <cstddef>
#include <concepts>
#include <span>

namespace uv::math::pde
{
	/**
	 * @brief Construct the mid-point (face) grid from a cell-centered grid.
	 *
	 * Returns a vector of size n-1 where each entry is the average of two
	 * consecutive points in xGrid. The input grid must be sorted.
	 */
	template <std::floating_point T>
	Vector<T> createXMidGrid(std::span<const T> xGrid);

	/**
	 * @brief Construct a discrete initial density for Fokker–Planck evolution.
	 *
	 * This function creates a probability density vector `p0` on a spatial grid
	 * suitable for forward (Fokker–Planck) PDE solvers.  The returned density
	 * approximates a Dirac delta located at the initial spatial point `x0`,
	 * ensuring that the total probability mass is approximately one.
	 *
	 *
	 * ### Behavior
	 * - Performs a binary search (`std::lower_bound`) to locate the closest grid
	 *   index to `x0`
	 * - Accounts for non-uniform grid spacing by estimating a local step size
	 * - Ensures total mass conservation under the discrete inner product
	 *
	 * @param x0     Initial spatial location where the density is concentrated
	 * @param xGrid  Sorted spatial grid defining the PDE domain
	 *
	 * @return A vector `p0` of the same size as `xGrid`, containing a discrete
	 *         probability density suitable for use as initial data in a
	 *         Fokker–Planck PDE solver.
	 */
	template <std::floating_point T>
	Vector<T> fokkerPlankInit(T x0,
	std::span<const T> xGrid);

	/**
	 * @brief Compute the Chang-Cooper interpolation weight at a cell interface.
	 *
	 * This function evaluates the Chang-Cooper weight
	 *
	 *   delta(w) = 1 / w - 1 / (exp(w) - 1)
	 *
	 * where
	 *
	 *   w = (B / C) * dx
	 *
	 * is the dimensionless drift-diffusion ratio at the cell interface.
	 *
	 * The weight delta in [0, 1] determines the upwind bias in the conservative
	 * finite-difference discretisation of the Fokker-Planck flux.
	 *
	 * Properties of the Chang-Cooper scheme:
	 *  - preserves non-negativity of the solution,
	 *  - conserves total probability mass,
	 *  - exactly reproduces the stationary (equilibrium) solution,
	 *  - smoothly interpolates between central differencing (w -> 0, delta -> 0.5)
	 *    and full upwinding (|w| -> infinity).
	 *
	 * Numerical safeguards:
	 *  - Non-finite values (NaN or Inf) are detected and replaced by delta = 0.5.
	 *  - The result is clamped to the theoretical bounds [0, 1].
	 *  - A warning is emitted if clamping occurs.
	 *
	 * Reference:
	 *   J. S. Chang and G. Cooper,
	 *   "A Practical Difference Scheme for Fokker-Planck Equations",
	 *   Journal of Computational Physics, Vol. 6, pp. 1-16, 1970.
	 *
	 * @param i Grid index in the first dimension.
	 * @param j Grid index in the second dimension.
	 * @param x Dimensionless drift-diffusion ratio w at the cell interface.
	 * @return Chang-Cooper weight delta in the interval [0, 1].
	 */
	template <std::floating_point T>
	core::Matrix<T> changCooperWeights(const core::Matrix<T>& B,
		const core::Matrix<T>& C,
        T dx);

	//core::Matrix<Real> fokkerPlankSolve(const Vector<Real>& pdeInitCond,
	//	const Vector<Real>& dTGrid,
	//	const Vector<Real>& dXGrid,
	//	const TriDiag& coefficients);

} // namespace uv::math::pde

#include "Functions.inl"