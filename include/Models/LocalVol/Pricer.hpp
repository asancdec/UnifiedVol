// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.hpp
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
#include "Core/MarketData.hpp"

#include <cstddef>
#include <array>

namespace uv::models::localvol
{
	class Pricer
	{
	private:

		// Data to price
		Vector<Real> tenors_;				// Tenors 
		Vector<Real> strikes_;				// Strikes

		// Market data
		Real S_;							// Spot price
		Real r_;							// Risk-free rate
		Real q_;							// Dividend yield

		// Grid dimensions
		const std::size_t NT_;				// Time points
		const std::size_t NS_;				// Spot points

		// PDE grids
		Vector<Real> timeGrid_;				// Time grid
		Vector<Real> spotGrid_;				// Spot grid


	public:

		// ---------- Constructors ----------

		Pricer() = delete;


		//explicit Pricer(const Vector<Real>& tenors,
		//	const Vector<Real>& strikes,
		//	const Vector<Real>& timeGrid_,
		//	const Vector<Real>& spotGrid_
		//);

		explicit Pricer(const Vector<Real>& tenors,
			const Vector<Real>& strikes,
			const uv::core::MarketData marketData,
			const std::size_t NT,
			const std::size_t NS,
			const unsigned int X = 3.0);

		// ---------- Pricing ----------

		Vector<Real> priceCall(const Matrix<Real>& localVol,
			const Vector<Real>& surfaceTenors,
			const Vector<Real>& surfaceStrikes);

	};
}  // uv::models::localvol