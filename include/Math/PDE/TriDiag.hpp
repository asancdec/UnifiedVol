// SPDX-License-Identifier: Apache-2.0
/*
 * File:        TriDiag.hpp
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
#include "Core/Matrix.hpp"

namespace uv::math::pde
{
	/**
	 * @brief Tridiagonal coefficient matrices for 1D PDE discretization.
	 *
	 * This struct stores three matrices corresponding to the sub-diagonal,
	 * main diagonal, and super-diagonal entries that arise in finite-difference
	 * discretization of 1D parabolic PDEs (e.g., Fokker–Planck or Black–Scholes).
	 *
	 * Each matrix has dimensions `NT x NS` where:
	 *   - `NT` is the number of time steps,
	 *   - `NS` is the number of spatial grid points.
	 *
	 * - `lower(i,j)` contains the coefficient multiplying the value at index `(i,j-1)`
	 * - `diag(i,j)`  contains the coefficient multiplying the value at index `(i,j)`
	 * - `upper(i,j)` contains the coefficient multiplying the value at index `(i,j+1)`
	 *
	 * These coefficients are used to assemble a tridiagonal linear system
	 * for implicit or Crank–Nicolson time stepping schemes.
	 */
	struct TriDiag
	{
		core::Matrix<Real> lower; ///< Sub-diagonal coefficients (A)
		core::Matrix<Real> diag;  ///< Main diagonal coefficients (B)
		core::Matrix<Real> upper; ///< Super-diagonal coefficients (C)
	};
}  // namespace uv::math::pde
