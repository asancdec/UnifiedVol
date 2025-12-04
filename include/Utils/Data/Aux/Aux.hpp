/*
* Aux.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <vector>

namespace uv
{
	// Transposes matrix (vector of vectors)
	template <typename T>
	std::vector<std::vector<T>> transposeMatrix(const std::vector<std::vector<T>>& input);
}

#include "Aux.inl"