// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
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

#include "Models/LocalVol/Pricer.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Math/Interpolation.hpp"
#include "Math/PDE/Functions.hpp"
#include "Core/Functions.hpp"

#include <algorithm>   
#include <string>    
#include <utility>
#include <cmath>

#include <iostream>

namespace uv::models::localvol
{
    // Erase later
    using namespace uv;
    using namespace uv::core;

    Pricer::Pricer(const Vector<Real>& tenors,
        const Vector<Real>& strikes,
        const uv::core::MarketData marketData,
        const std::size_t NT,
        const std::size_t NS,
        const unsigned int X)
        : tenors_(tenors),
        strikes_(strikes),
        S_(marketData.S),
        r_(marketData.r),
        q_(marketData.q),
        NT_(NT),
        NS_(NS)
    {
        validateInputs();
        initPDEMainGrids(X);
        initPDECachedGrids();
    }

    Vector<Real> Pricer::price(const Matrix<Real>& localVar,
        const Vector<Real>& surfaceTenors,
        const Vector<Real>& surfaceStrikes)
    {
        // ---------- Interpolate local vol ----------

        // Tensor spline PCHIP product
        const Matrix<Real> lvGrid
        {
            uv::math::interp::pchipInterp2D
            (
                spotGrid_,
                timeGrid_,
                surfaceStrikes,
                surfaceTenors,
                localVar
            )
        };

        // ---------- Fokker Plank PDE ----------

        // Complete diffusion term with local var
        Matrix<Real> diffusion
        {
            uv::core::hadamard
            (
                lvGrid,
                pdeDiffusion_
            )
        };

        // Solve the PDE
        Matrix<Real> pdf
        {
            uv::math::pde::fokkerPlankSolve
            (
               pdeInitCond_,
               dTGrid_,
               dSGrid_,
               pdeDrift_,
               diffusion
            )
        };
        

        return surfaceTenors;
    }

    void Pricer::validateInputs()
    {
       // Sorted tenors
        UV_REQUIRE(
            std::is_sorted(tenors_.begin(), tenors_.end()),
            ErrorCode::InvalidArgument,
            "Pricer: tenors must be non-decreasing"
        );

        // Sorted strikes
        UV_REQUIRE(
            std::is_sorted(strikes_.begin(), strikes_.end()),
            ErrorCode::InvalidArgument,
            "Pricer: strikes must be strictly increasing"
        );
    }

    void Pricer::initPDEMainGrids(const Real X)
    {
        // Time axis
        const Real maxT{ tenors_.back() };
        timeGrid_ = uv::core::generateGrid(Real(0.0), maxT, NT_);
        dTGrid_ = uv::core::diff(timeGrid_);

        // Forward axis
        const Real maxS{ S_ * Real(X) }; 
        spotGrid_ = uv::core::generateGrid(Real(0.0), maxS, NS_);
        dSGrid_ = uv::core::diff(spotGrid_);
    }

    void Pricer::initPDECachedGrids()
    {
        // Initial condition
        // Approximation to pdf at t=0
        pdeInitCond_ = uv::math::pde::fokkerPlankInit(S_, spotGrid_);

        // Drift: (r - q) * S_i
        pdeDrift_ = Matrix<Real>
        (
            NT_,
            uv::core::multiply(spotGrid_, (r_ - q_))
        );

        // Diffusion: 0.5 * S_i^2
        pdeDiffusion_ = uv::core::multiply
        (
            uv::core::hadamard(spotGrid_, spotGrid_),
            Real(0.5)
        );
    }




} // namespace uv::models::localvol