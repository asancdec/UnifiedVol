// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.cpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include "Utils/Aux/Errors.hpp"  
#include "Core/Functions.hpp"

#include <string> 

namespace uv::core
{
    Matrix<Real> transposeMatrix(const Matrix<Real>& input)
	{
        // ---------- Check matching dimensions ----------

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

        // ---------- Transpose matrix ----------
        
        std::vector<Vector<Real>> result(numCols, Vector<Real>(numRows));

        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                result[j][i] = input[i][j];
            }
        }

        return result;
	}

    Vector<Real> generateGrid(const Real bound1,
        const Real bound2,
        const size_t steps
    ) noexcept
    {
        // ---------- Allocate memory ----------

        Vector<Real> grid;
        grid.reserve(steps);

        // ---------- Handle special case ----------

        if (steps <= 1)
        {
            grid.push_back(static_cast<Real>(bound1));
            return grid;
        }

        // ---------- Generate grid ----------

        const Real dx{ (bound2 - bound1) / static_cast<Real>(steps - 1) };

        for (std::size_t i = 0; i < steps; ++i)
        {
            grid.push_back(static_cast<Real>(bound1) + dx * static_cast<Real>(i));
        }

        return grid;
    }

} // namespace uv::core