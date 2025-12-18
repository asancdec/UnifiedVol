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

#include "Utils/Types.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Math/PDE/TriDiag.hpp"
#include "Core/Matrix/Matrix.hpp"

#include <cstddef>

namespace uv::math::pde
{
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
	Vector<Real> fokkerPlankInit(Real x0,
	const Vector<Real>& xGrid);

	core::Matrix<Real> fokkerPlankSolve(const Vector<Real>& pdeInitCond,
		const Vector<Real>& dTGrid,
		const Vector<Real>& dXGrid,
		const TriDiag& coefficients);

} // namespace uv::math::pde