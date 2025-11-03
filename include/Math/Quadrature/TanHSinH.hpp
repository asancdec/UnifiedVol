/**
* TanHSinH.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <array>
#include <cstddef>

namespace uv
{
	template <::std::size_t N>
	class TanHSinH
	{
	private:

		//--------------------------------------------------------------------------
		// Internal struct
		//--------------------------------------------------------------------------	
		struct Node
		{
			long double y;     // yn := 1 - xn
			long double x;     // Abscissas value
			long double w;     // Weight value
		};

		//--------------------------------------------------------------------------
		// Member variables
		//--------------------------------------------------------------------------	

		const long double     h_;            // Step size
		::std::array<Node, N> nodes_;        // Nodes vector

		//--------------------------------------------------------------------------
		// Math
		//--------------------------------------------------------------------------

		// Fixed Tanh-Sinh rule node (returns struct with weight and abscissas value)
		static Node generateNode(const long double nh) noexcept;

		// Domain mapping function (0, inf) -> (-1, 1)
		template<typename F>
		static long double transformIntegrand(long double x, long double y, F&& f) noexcept;

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

		//--------------------------------------------------------------------------
		// Utilities
		//--------------------------------------------------------------------------
		void printGrid() const noexcept;
	};
}

#include "TanHSinH.inl"
