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
#include "Core/Functions.hpp"

#include <algorithm>   
#include <string>     

#include <iostream>

namespace uv::models::localvol
{
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
        // ---------- Check sorted tenors ----------

        UV_REQUIRE(
            std::is_sorted(tenors_.begin(), tenors_.end()),
            ErrorCode::InvalidArgument,
            "Pricer: tenors must be non-decreasing"
        );

        // ---------- Check sorted strikes ----------
        
        UV_REQUIRE(
            std::is_sorted(strikes_.begin(), strikes_.end()),
            ErrorCode::InvalidArgument,
            "Pricer: strikes must be strictly increasing"
        );

        // ---------- Generate grids ----------

        timeGrid_ = uv::core::generateGrid(Real(0.0), tenors.back(), NT);
        spotGrid_ = uv::core::generateGrid(Real(0.0), S_ * Real(X), NS);
    }

    Vector<Real> Pricer::priceCall(const Matrix<Real>& localVol,
        const Vector<Real>& surfaceTenors,
        const Vector<Real>& surfaceStrikes)
    {
        // ---------- Interpolate local vol ----------

        const Matrix<Real> lvGrid{};
        
        return surfaceTenors;
    }


} // namespace uv::models::localvol