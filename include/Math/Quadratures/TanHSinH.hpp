// SPDX-License-Identifier: Apache-2.0
/*
 * File:        TanHSinH.hpp
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

#include <array>
#include <cstddef>

namespace uv::math
{
	template <std::size_t N>
	class TanHSinH
	{
	private:

		//--------------------------------------------------------------------------
		// Internal struct
		//--------------------------------------------------------------------------	
		struct Node
		{
			Real w;              // weight value
			Real y;              // yn term
			Real x;              // abscissas value
			Real factorRight;    // scaling term RHS
			Real inputRight;     // transformed input RHS
			Real factorLeft;     // scaling term LHS
			Real inputLeft;      // transformed input LHS
		};

      	//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------	

		const Real       h_;            // Step size
		std::array<Node, N> nodes_;     // Nodes vector

		//--------------------------------------------------------------------------
		// Math
		//--------------------------------------------------------------------------

		// Fixed Tanh-Sinh rule node (returns struct with weight and abscissas value)
		Node generateNode(const Real nh) const noexcept;

	public:

		//--------------------------------------------------------------------------
		// Initialization
		//--------------------------------------------------------------------------
		TanHSinH();

		//--------------------------------------------------------------------------
		// Math
		//--------------------------------------------------------------------------

		// Evaluation a function's numeric integral
		template<typename F>
		Real integrateZeroToInf(F&& f) const noexcept;

		// Evaluation of multiple numeric integrals
		template<std::size_t M, typename F >
		std::array<Real, M> integrateZeroToInfMulti(F&& f) const noexcept;

		//--------------------------------------------------------------------------
		// Utilities
		//--------------------------------------------------------------------------
		void printGrid() const noexcept;
	};
}

#include "TanHSinH.inl"
