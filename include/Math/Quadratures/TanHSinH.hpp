/**
* TanHSinH.hpp
* Author: Alvaro Sanchez de Carlos
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
