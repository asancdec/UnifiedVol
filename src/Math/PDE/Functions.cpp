// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
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

#include "Math/PDE/Functions.hpp"
#include "Utils/Aux/Errors.hpp"

#include <cstddef>
#include <cmath>
#include <string>
#include <algorithm>

namespace uv::math::pde
{
    Vector<Real> fokkerPlankInit(const Real x0,
        const Vector<Real>& xGrid)
    {	
        // ---------- Check size ----------

		std::size_t NX{ xGrid.size() };

        // Throw if grid size is one or less
        UV_REQUIRE(
            NX > 1,
            ErrorCode::InvalidArgument,
            "fokkerPlankInit: grid size must be > 1 (got " + std::to_string(NX) + ")"
        );

        // ---------- Check sorted ----------

        // Throw if not sorted
        UV_REQUIRE(
            std::is_sorted(xGrid.begin(), xGrid.end()),
            ErrorCode::InvalidArgument,
            "fokkerPlankInit: xGrid must be non-decreasing"
        );

        // ---------- Check in bounds ----------

        // Throw if out of bounds
        UV_REQUIRE(
            x0 >= xGrid.front() && x0 <= xGrid.back(),
            ErrorCode::InvalidArgument,
            "fokkerPlankInit: Initial value must lie within the grid domain"
        );

        // ---------- Binary search index ----------

        auto it = std::lower_bound(xGrid.begin(), xGrid.end(), x0);
        std::size_t idx{ static_cast<std::size_t>(std::distance(xGrid.begin(), it)) };

        // If there is an element to the left
        if (idx > 0)
        {
            // Choose the closest of idx and idx-1
            if (std::abs(xGrid[idx] - x0) > std::abs(xGrid[idx - 1] - x0))
            {
                --idx;
            }
        }
        
        // ---------- Dirac approximation ----------

        Vector<Real> p0(NX, Real(0.0));
        Real dX{};

        // Step size at the left 
        if (idx == 0)
        {
            dX = xGrid[1] - xGrid[0];
        }
        // Step size at the right
        else if (idx == NX - 1)
        {
            dX = xGrid[NX - 1] - xGrid[NX - 2];
        }
        else
        {   
            // Average step sizes
            const Real dXLeft{ xGrid[idx] - xGrid[idx - 1] };
            const Real dXRight{ xGrid[idx + 1] - xGrid[idx] };
            dX = Real(0.5) * (dXLeft + dXRight);
        }

        // Probability mass under p0 is 1
        // Therefore the following must hold approximately
        p0[idx] = Real(1.0) / dX;

        return p0;
    }

    //Matrix<Real> fokkerPlankSolve(const Vector<Real>& pdeInitCond,
    //    const Vector<Real>& dTGrid,
    //    const Vector<Real>& dXGrid,
    //    const TriDiag& coefficients)
    //{
    //    // ---------- Dimension checks ----------

    //    const std::size_t nT{ drift.size() };
    //    const std::size_t nX{ drift[0].size() };

    //    UV_REQUIRE(
    //        diffusion.size() == nT,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion must have same number of rows as drift"
    //    );

    //    UV_REQUIRE(
    //        diffusion[0].size() == nX,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion must have same number of columns as drift"
    //    );

    //    UV_REQUIRE(
    //        pdeInitCond.size() == nX,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: pdeInitCond size must equal number of spatial points (got " +
    //        std::to_string(pdeInitCond.size()) + " and " + std::to_string(nX) + ")"
    //    );

    //    UV_REQUIRE(
    //        dTGrid.size() == (nT - 1),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: dTGrid size must be nTimeLevels - 1 (got " +
    //        std::to_string(dTGrid.size()) + " and " + std::to_string(nT - 1) + ")"
    //    );

    //    UV_REQUIRE(
    //        dXGrid.size() == (nX - 1),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: dXGrid size must be nSpacePoints - 1 (got " +
    //        std::to_string(dXGrid.size()) + " and " + std::to_string(nX - 1) + ")"
    //    );

    //    // ---------- Math checks ----------

    //    UV_REQUIRE(
    //        std::all_of(dTGrid.begin(), dTGrid.end(),
    //            [](Real dt) { return dt > Real(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all time steps dTGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(dXGrid.begin(), dXGrid.end(),
    //            [](Real dx) { return dx > Real(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all space steps dXGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(diffusion.begin(), diffusion.end(),
    //            [](const Vector<Real>& row)
    //            {
    //                return std::all_of(row.begin(), row.end(),
    //                    [](Real v) { return v >= Real(0); });
    //            }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion coefficients must be non-negative"
    //    );

    //    // ---------- Internal coefficients ----------




    //    return diffusion;
    //}
    //Matrix<Real> fokkerPlankSolve(const Vector<Real>& pdeInitCond,
    //    const Vector<Real>& dTGrid,
    //    const Vector<Real>& dXGrid,
    //    const TriDiag& coefficients)
    //{
    //    // ---------- Dimension checks ----------

    //    const std::size_t nT{ drift.size() };
    //    const std::size_t nX{ drift[0].size() };

    //    UV_REQUIRE(
    //        diffusion.size() == nT,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion must have same number of rows as drift"
    //    );

    //    UV_REQUIRE(
    //        diffusion[0].size() == nX,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion must have same number of columns as drift"
    //    );

    //    UV_REQUIRE(
    //        pdeInitCond.size() == nX,
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: pdeInitCond size must equal number of spatial points (got " +
    //        std::to_string(pdeInitCond.size()) + " and " + std::to_string(nX) + ")"
    //    );

    //    UV_REQUIRE(
    //        dTGrid.size() == (nT - 1),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: dTGrid size must be nTimeLevels - 1 (got " +
    //        std::to_string(dTGrid.size()) + " and " + std::to_string(nT - 1) + ")"
    //    );

    //    UV_REQUIRE(
    //        dXGrid.size() == (nX - 1),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: dXGrid size must be nSpacePoints - 1 (got " +
    //        std::to_string(dXGrid.size()) + " and " + std::to_string(nX - 1) + ")"
    //    );

    //    // ---------- Math checks ----------

    //    UV_REQUIRE(
    //        std::all_of(dTGrid.begin(), dTGrid.end(),
    //            [](Real dt) { return dt > Real(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all time steps dTGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(dXGrid.begin(), dXGrid.end(),
    //            [](Real dx) { return dx > Real(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all space steps dXGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(diffusion.begin(), diffusion.end(),
    //            [](const Vector<Real>& row)
    //            {
    //                return std::all_of(row.begin(), row.end(),
    //                    [](Real v) { return v >= Real(0); });
    //            }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion coefficients must be non-negative"
    //    );

    //    // ---------- Internal coefficients ----------




    //    return diffusion;
    //}


} // namespace uv::math::pde