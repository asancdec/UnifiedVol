/*
* Helpers.hpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Aux/Errors.hpp"    

#include <string>

namespace uv::utils
{
	template <typename T>
	std::vector<std::vector<T>> transposeMatrix(const std::vector<std::vector<T>>& input)
	{
        // ---------- Sanity checks ------

        // Extract dimensions
        const std::size_t numRows{ input.size() };
        const std::size_t numCols{ input[0].size() };

        // Validate that all rows have the same number of columns
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
        
        // Allocate result with swapped dimensions
        // result[j][i] = input[i][j]
        std::vector<std::vector<T>> transposed(numCols, std::vector<T>(numRows));

        // Perform transpose
        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                transposed[j][i] = input[i][j];
            }
        }

        return transposed;
	}
}