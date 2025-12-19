// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.cpp
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


#include "Models/SVI/Functions.hpp"  
#include "Core/VolSurface.hpp"
#include "Core/Matrix/Functions.hpp"
#include "Models/SVI/Params.hpp"
#include "Utils/Aux/Errors.hpp"

#include <array>
#include <cmath>      
#include <cstddef>
#include <span>


namespace uv::models::svi::detail
{
    std::array<double, 4> initGuess() noexcept
    {
        return 
        { 
            0.1,        // b
            -0.5,       // rho
            0.1,        // m
            0.1         // sigma
        };
    }

    std::array<double, 4> lowerBounds(double logKFMin) noexcept
    {
        return
        {
            0.001,                               // b            
            -0.9999,                             // rho
            5.0 * logKFMin,                      // m
            0.01                                 // sigma
        };
    }

    std::array<double, 4> upperBounds(double logKFMax) noexcept
    {
        return
        {
            2.0,                                // b
            0.9999,                             // rho
            5.0 * logKFMax,                     // m
            10.0                                // sigma
        };
    }

    double aParam(double atmWK,
        double b,
        double rho,
        double m,
        double sigma) noexcept
    {
        return atmWK - b * (-rho * m + std::sqrt(m * m + sigma * sigma));
    }

    double gk(const GkCache& p) noexcept
    {
        // g(k) = A^2 - B * (w')^2 / 4 + w''/2
        return (p.A * p.A) - 0.25 * p.wkD1Squared * p.B + p.wkD2 * 0.5;  
    }

    std::array<double, 4> gkGrad(
        double b,
        double rho,
        double m,
        double sigma,
        double k,
        const GkCache& p
    ) noexcept
    {
        // ---------- Extract variables ----------

        const double x{ p.x };
        const double R{ p.R };
        const double invR{ p.invR };
        const double invRCubed{ p.invRCubed };
        const double invR5{ p.invR5 };
        const double sigmaSquared{ p.sigmaSquared };
        const double wkInv{ p.wkInv };
        const double wkSquaredInv{ p.wkSquaredInv };
        const double wkD1{ p.wkD1 };
        const double wkD1Squared{ p.wkD1Squared };
        const double A{ p.A };
        const double B{ p.B };
        const double R0{ p.R0 };
        const double invR0{ p.invR0 };

        // ---------- ∂g/∂w, ∂g/∂w', ∂g/∂w'' ----------

        const double dgdw{ A * k * wkD1 * wkSquaredInv + 0.25 * wkD1Squared * wkSquaredInv };
        const double dgdw1{ -A * k * wkInv - 0.5 * wkD1 * B };
        const double dgdw2{ 0.5 };

        // ---------- ∂w/∂θ ----------

        std::array<double, 4> dw
        {
            (rho * x + R) + (rho * m - R0),                   
            (b * x) + (b * m),                                 
            (-b * (rho + x * invR)) + (b * (rho - m * invR0)), 
            (b * sigma * invR) + (-b * (sigma * invR0))      
        };

        // ---------- ∂w′/∂θ ----------

        std::array<double, 4> dw1
        {
            rho + x * invR,
            b,
            -b * sigmaSquared * invRCubed,
            -b * x * sigma * invRCubed
        };

        // ---------- ∂w″/∂θ ----------

        std::array<double, 4> dw2
        {
            sigmaSquared * invRCubed,
            0.0,
            3.0 * b * sigmaSquared * x * invR5,
            b * (2.0 * sigma * invRCubed - 3.0 * sigmaSquared * sigma * invR5)
        };

        // ---------- Chain rule ----------

        std::array<double, 4> dg{};
        for (int j = 0; j < 4; ++j)
        {
            dg[j] = dgdw * dw[j] + dgdw1 * dw1[j] + dgdw2 * dw2[j];
        }
        return dg;
    }
} // namespace uv::models::svi::detail