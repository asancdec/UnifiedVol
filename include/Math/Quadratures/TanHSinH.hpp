/**
* TanHSinH.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

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
			long double w;              // weight value
			long double y;              // yn term
			long double x;              // abscissas value
			long double factorRight;    // scaling term RHS
			long double inputRight;     // transformed input RHS
			long double factorLeft;     // scaling term LHS
			long double inputLeft;      // transformed input LHS
		};

      	//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------	

		const long double     h_;            // Step size
		std::array<Node, N> nodes_;        // Nodes vector

		//--------------------------------------------------------------------------
		// Math
		//--------------------------------------------------------------------------

		// Fixed Tanh-Sinh rule node (returns struct with weight and abscissas value)
		Node generateNode(const long double nh) const noexcept;

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
		long double integrateZeroToInf(F&& f) const noexcept;

		// Evaluation of multiple numeric integrals
		template<std::size_t M, typename F >
		std::array<long double, M> integrateZeroToInfMulti(F&& f) const noexcept;

		//--------------------------------------------------------------------------
		// Utilities
		//--------------------------------------------------------------------------
		void printGrid() const noexcept;
	};
}

#include "TanHSinH.inl"
