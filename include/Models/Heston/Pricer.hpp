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
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Config.hpp"
#include "Models/Heston/CharFunData.hpp"
#include "Math/Quadratures/TanHSinH.hpp"
#include "Core/VolSurface.hpp"

#include <memory>
#include <array>
#include <optional>    
#include <cstddef>  
#include <tuple>

namespace uv::models::heston
{	

	template <std::size_t N>
	class Pricer
	{
	private:

		//--------------------------------------------------------------------------
		// Private member variables
		//--------------------------------------------------------------------------
		std::optional<Params> params_;
		std::shared_ptr<const math::TanHSinH<N>> quad_;
		const Config config_;

		//--------------------------------------------------------------------------
		// Pricing
		//--------------------------------------------------------------------------
		// Calculate residues arising from the contour shift
		static Real getResidues(Real alpha,
			Real F,
			Real K) noexcept;

		// Get ITM or OTM damping parameter
		Real getAlpha(Real w) const noexcept;

		// Determine contour shift angle
		static Real getPhi(Real kappa,
			Real theta,
			Real sigma,
			Real rho,
			Real v0,
			Real T,
			Real w) noexcept;

		// Returns only price
		// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
		static Complex<Real> charFunction(Real kappa,
			Real theta,
			Real sigma,
			Real rho,
			Real v0,
			Real T,
			const Complex<Real>& u) noexcept;

		//--------------------------------------------------------------------------
		// Calibration
		//--------------------------------------------------------------------------
		// Returns a struct of precalculated variables for efficient gradient computation
		// Albrecher, H., P. Mayer, W. Schoutens, and J. Tistaert (2007)
		static CharFunData charFunctionCal(Real  kappa,
			Real theta,
			Real sigma,
			Real rho,
			Real v0,
			Real T,
			const Complex<Real>& u) noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		Pricer() = delete;
		explicit Pricer(std::shared_ptr<const math::TanHSinH<N>> quad,
			const Config& config);

		//--------------------------------------------------------------------------
		// Pricing
		//--------------------------------------------------------------------------
		// Overload 1: using user-defined parameters
		// Andersen & Lake Implementation
		Real callPrice(Real kappa,
			Real theta,
			Real sigma,
			Real rho,
			Real v0,
			Real T,
			Real F,
			Real r,
			Real K) const noexcept;

		// Overload 2: using class instance parameters
		// Andersen & Lake Implementation
		Real callPrice(Real T,
			Real F,
			Real r,
			Real K) const;

		// Calculate price and parameter gradient for the calibration
		std::array<Real, 6> callPriceWithGradient(Real kappa,
			Real theta,
			Real sigma,
			Real rho,
			Real v0,
			Real T,
			Real F,
			Real r,
			Real K) const noexcept;

		//--------------------------------------------------------------------------
		// Setters
		//--------------------------------------------------------------------------
		template <typename T>
		void setParams(const std::array<T, 5>& params) noexcept;
	};
} // namespace uv::models::heston

#include "Pricer.inl"