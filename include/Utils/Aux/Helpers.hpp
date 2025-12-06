/*
* Helpers.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

namespace uv::utils
{
    /**
     * @brief Transpose a dense matrix stored as a vector of row vectors.
     *
     * @tparam T        Element type.
     * @param  input    Matrix with shape [numRows][numCols].
     * @return          Transposed matrix with shape [numCols][numRows].
     */
	template <typename T>
	Matrix<T> transposeMatrix(const Matrix<T>& input);

}

#include "Helpers.inl"