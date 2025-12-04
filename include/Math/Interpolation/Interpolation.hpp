/**
* Interpolation.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <vector>

namespace uv
{
	// Returns forward differences (segment slopes) dy/dx on each interval
	// NOTE: if "extrapolateEnd" is true, the code will flat extrapolate the lest element
    template <typename T>
	std::vector<T> d1PieceWise(
		const std::vector<T>& x,
		const std::vector<T>& y,
		bool extrapolateEnd = true
	);
}

#include "Interpolation.inl"