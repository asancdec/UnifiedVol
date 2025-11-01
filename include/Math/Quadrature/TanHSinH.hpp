/**
* TanHSinH.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once 
#include <vector>
#include <cmath> 
#include <utility> 

namespace uv
{
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

		const long double h_;            // Step size
		::std::vector<Node> nodes_;        // Nodes vector
		::std::size_t       N_;            // Number of nodes


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
		TanHSinH(const unsigned int N = 64);

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
