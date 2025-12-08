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
		// Internal struct to save intermediate calculations
		//--------------------------------------------------------------------------
		struct CharFunCache
		{
			// Original CF outputs
			Complex<Real> psi;                 // psi(u) := exp( A + v0 * B )
			Complex<Real> A;                   // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
			Complex<Real> B;                   // B := u(u+i)*y / (1 − r*y)
			Complex<Real> beta;                // beta := kappa − sigma*rho*(i*u)
			Complex<Real> D;                   // D := sqrt( beta^2 + sigma^2*u(u+i) )
			Complex<Real> DT;                  // DT := D*T
			Complex<Real> betaPlusD;           // beta + D
			Complex<Real> betaMinusD;          // beta − D  (stable)
			Complex<Real> ui;                  // ui := u*i
			Complex<Real> kFac;                // kFac := (kappa*theta)/sigma^2
			Real invSigma2;					   // 1 / sigma^2
			Real kappaTheta;				   // kappa * theta
			Real sigma2;					   // sigma^2

			// Rescued intermediates (for gradient path)
			Complex<Real> uu;                  // uu := u(u+i)
			Complex<Real> eDT;                 // eDT := exp( −DT )
			Complex<Real> g;                   // g := (beta−D)/(beta+D)
			Complex<Real> Q;                   // Q := 1 − g*eDT
			Complex<Real> invQ;                // 1 / Q
			Complex<Real> invQ2;               // 1 / Q^2
			Complex<Real> R;                   // R := 1 − g
			Complex<Real> S;                   // S := (beta−D)*T − 2*log(Q/R)
			Complex<Real> fracB;               // fracB := (1 − eDT) / Q
			Complex<Real> denomG;              // denomG := (beta + D)^2
			Complex<Real> betaMinusDinvSigma2; // (betaMinusD)/sigma^2
		};

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
		static CharFunCache charFunctionCal(Real  kappa,
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