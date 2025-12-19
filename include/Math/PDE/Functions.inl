// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
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

#include "Core/Matrix/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <format>

namespace uv::math::pde
{

    template <std::floating_point T>
    Vector<T> createXMidGrid(std::span<const T> xGrid)
    {
        // ---------- Check size ----------

        const std::size_t n{ xGrid.size() };

        // Throw if less than two elements
        UV_REQUIRE
        (
            n >= 2,
            ErrorCode::InvalidArgument,
            "createXMidGrid: grid must have at least 2 points"
        );


        // ---------- Check sorted ----------

        // Throw if not sorted
        UV_REQUIRE(
            std::is_sorted(xGrid.begin(), xGrid.end()),
            ErrorCode::InvalidArgument,
            "createXMidGrid: xGrid must be non-decreasing"
        );

        // ---------- Allocate ----------

        Vector<T> xMidGrid(n - 1);

        for (std::size_t i = 0; i < n - 1; ++i)
        {
            xMidGrid[i] = T{ 0.5 } *(xGrid[i] + xGrid[i + 1]);
        }

        return xMidGrid;
    }

    template <std::floating_point T>
    Vector<T> fokkerPlankInit(const T x0,
        std::span<const T> xGrid)
    {
        // ---------- Check size ----------

        std::size_t nx{ xGrid.size() };

        // Throw if grid size is one or less
        UV_REQUIRE(
            nx > 1,
            ErrorCode::InvalidArgument,
            "fokkerPlankInit: grid size must be > 1 (got " + std::to_string(nx) + ")"
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

        Vector<T> p0(nx, T(0));
        T dX;

        // Step size at the left 
        if (idx == 0)
        {
            dX = xGrid[1] - xGrid[0];
        }
        // Step size at the right
        else if (idx == nx - 1)
        {
            dX = xGrid[nx - 1] - xGrid[nx - 2];
        }
        else
        {
            // Average step sizes
            const T dXLeft{ xGrid[idx] - xGrid[idx - 1] };
            const T dXRight{ xGrid[idx + 1] - xGrid[idx] };
            dX = T(0.5) * (dXLeft + dXRight);
        }

        // Probability mass under p0 is 1
        // Therefore the following must hold approximately
        p0[idx] = T(1) / dX;

        return p0;
    }

    template <std::floating_point T>
    core::Matrix<T> changCooperWeights(const core::Matrix<T>& B,
        const core::Matrix<T>& C,
        T dx)
    {
        // ---------- Dimension checks ----------

        UV_REQUIRE
        (
            B.rows() == C.rows() &&
            B.cols() == C.cols(),
            ErrorCode::InvalidArgument,
            "changCooperWeights: invalid dimensions"
        );

        // ---------- Calculate weights ----------

        core::Matrix<T> w{ core::divide(B, C) };
        w *= dx;

        core::Matrix<T> weights
        {
            core::transformIndexed<T>
            (
                w,
               [](std::size_t i, std::size_t j, T x)
                {
                    // ---------- Chang-Cooper weight ----------

                    // expm1(x) = exp(x)-1
                    T delta {T{1} / x - T{1} / std::expm1(x)};

                    // ---------- Handle NaN/Inf ----------

                    if (!std::isfinite(delta))
                    {
                        UV_WARN(
                            true,
                            std::format(
                                "ChangCooperWeights: non-finite delta at ({}, {}), "
                                "omega = {:.6e} -> forcing 0.5",
                                i, j, x
                            )
                        );
                        return T{ 0.5 };
                    }

                    const T clamped{ std::clamp(delta, T{ 0 }, T{ 1 }) };

                    // Warn only if clamping occurred
                    UV_WARN(
                        clamped != delta,
                        std::format(
                            "ChangCooperWeights: delta clamped at ({}, {}): "
                            "omega = {:.6e}, delta = {:.6e} -> {:.6e}",
                            i, j, x, delta, clamped
                        )
                    );

                    return clamped;
                }
            )
        };

        return weights;
    }

    //Matrix<T> fokkerPlankSolve(const Vector<T>& pdeInitCond,
    //    const Vector<T>& dTGrid,
    //    const Vector<T>& dXGrid,
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
    //            [](T dt) { return dt > T(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all time steps dTGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(dXGrid.begin(), dXGrid.end(),
    //            [](T dx) { return dx > T(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all space steps dXGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(diffusion.begin(), diffusion.end(),
    //            [](const Vector<T>& row)
    //            {
    //                return std::all_of(row.begin(), row.end(),
    //                    [](T v) { return v >= T(0); });
    //            }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion coefficients must be non-negative"
    //    );

    //    // ---------- Internal coefficients ----------




    //    return diffusion;
    //}
    //Matrix<T> fokkerPlankSolve(const Vector<T>& pdeInitCond,
    //    const Vector<T>& dTGrid,
    //    const Vector<T>& dXGrid,
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
    //            [](T dt) { return dt > T(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all time steps dTGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(dXGrid.begin(), dXGrid.end(),
    //            [](T dx) { return dx > T(0); }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: all space steps dXGrid must be strictly positive"
    //    );

    //    UV_REQUIRE(
    //        std::all_of(diffusion.begin(), diffusion.end(),
    //            [](const Vector<T>& row)
    //            {
    //                return std::all_of(row.begin(), row.end(),
    //                    [](T v) { return v >= T(0); });
    //            }),
    //        ErrorCode::InvalidArgument,
    //        "fokkerPlankSolve: diffusion coefficients must be non-negative"
    //    );

    //    // ---------- Internal coefficients ----------




    //    return diffusion;
    //}


} // namespace uv::math::pde