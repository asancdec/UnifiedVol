/*
* Helpers.hpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Aux/Errors.hpp"   

#include <string> 

namespace uv::utils
{
	template <typename T>
    Matrix<T> transposeMatrix(const Matrix<T>& input)
	{
        // ---------- Check matching dimensions ------

        const std::size_t numRows{ input.size() };
        const std::size_t numCols{ input[0].size() };

        // Throw if all rows are not the same size
        for (std::size_t i = 1; i < numRows; ++i)
        {
            const std::size_t actual{ input[i].size() };
            const std::size_t expected{ numCols };

            UV_REQUIRE(
                actual == expected,
                ErrorCode::InvalidArgument,
                "transposeVector: inconsistent row length — row " +
                std::to_string(i) + " has " + std::to_string(actual) +
                " elements, expected " + std::to_string(expected)
            );
        }

        // ---------- Transpose matrix ------
        
        std::vector<Vector<T>> result(numCols, Vector<T>(numRows));

        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                result[j][i] = input[i][j];
            }
        }

        return result;
	}
} // namespace uv::utils