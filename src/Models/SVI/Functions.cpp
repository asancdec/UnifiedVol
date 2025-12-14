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
#include "Utils/Aux/Errors.hpp"

#include <cmath>      

namespace uv::models::svi
{
    core::VolSurface buildSurface(const core::VolSurface& volSurface,
        const Vector<Params>& params)
    {
        // ---------- Extract data ----------

        const std::size_t numTenors{ volSurface.numTenors() };
        const std::size_t numStrikes{ volSurface.numStrikes() };
        const Matrix<Real> kMatrix{ volSurface.logKFMatrix() };

        // ---------- Validate inputs ----------

        UV_REQUIRE(
            params.size() == numTenors,
            ErrorCode::InvalidArgument,
            "buildSurface: params size must equal number of tenors"
        );

        // ---------- Build total variance matrix ----------

        Matrix<Real> wK(numTenors, Vector<Real>(numStrikes));

        for (std::size_t i = 0; i < numTenors; ++i)
        {
            // Extract data
            const Params& p{ params[i] };
            const Vector<Real>& kSlice{kMatrix[i]};
            const double a{ double(p.a) };
            const double b{ double(p.b) };
            const double rho{ double(p.rho) };
            const double m{ double(p.m) };
            const double sigma{ double(p.sigma) };

            // Calculate total variance per strike
            for (std::size_t j = 0; j < numStrikes; ++j)
            {
                wK[i][j] = Real(detail::calculateWk(a, b, rho, m, sigma, kSlice[j]));
            }
        }

        // ---------- Build surface  ---------- 
        
        // Copy surface
        core::VolSurface sviVolSurface{ volSurface };

        // Set data
        sviVolSurface.setWt(wK);

        return sviVolSurface;
    }

    double gk(Real a, Real b, Real rho, Real m, Real sigma, Real k) noexcept
    {
        Real x{ k - m };                                          // x := k-m
        Real R{ std::hypot(x, sigma) };                           // R:= sqrt(x^2 + sigma^2)
        Real invR{ 1.0 / R };                                     // invR := 1 / R
        Real wk{ std::fma(b, (rho * x + R), a) };                 // w(k) = a + b*(rho*x + R)
        Real wkD1{ b * (rho + x * invR) };                        // w'(k) = b * (rho + x/R)
        Real wkD1Squared{ wkD1 * wkD1 };                          // w'(k)^2
        Real invRCubed{ invR * invR * invR };                     // 1/(R^3)
        Real sigmaSquared{ sigma * sigma };                       // sigma^2
        Real wkD2{ b * sigmaSquared * invRCubed };                // w''(k) = b * sigma^2 / R^3
        Real A{ 1.0 - 0.5 * k * wkD1 / wk };                      // A := 1 - k * w'/(2 * w)                        
        Real B{ (1.0 / wk) + 0.25 };                              // B := 1/w(k) + 1/4

        return (A * A) - 0.25 * wkD1Squared * B + wkD2 * 0.5;
    }
} // namespace uv::models::svi

namespace uv::models::svi::detail
{
    void validateInputs(const Vector<Real>& tenors,
        const Matrix<Real>& kMatrix,
        const Matrix<Real>& wMMatrix)
    {
        UV_REQUIRE(
            !tenors.empty(),
            ErrorCode::InvalidArgument,
            "validateInputs: tenors is empty"
        );

        UV_REQUIRE(
            !kMatrix.empty(),
            ErrorCode::InvalidArgument,
            "validateInputs: kMatrix is empty"
        );

        UV_REQUIRE(
            !wMMatrix.empty(),
            ErrorCode::InvalidArgument,
            "validateInputs: wMMatrix is empty"
        );

        UV_REQUIRE(
            kMatrix.size() == tenors.size(),
            ErrorCode::InvalidArgument,
            "validateInputs: kMatrix rows must equal number of tenors"
        );

        UV_REQUIRE(
            wMMatrix.size() == tenors.size(),
            ErrorCode::InvalidArgument,
            "validateInputs: wMMatrix rows must equal number of tenors"
        );

        UV_REQUIRE(
            !kMatrix[0].empty(),
            ErrorCode::InvalidArgument,
            "validateInputs: kMatrix has zero strikes"
        );

        const std::size_t numTenors{ tenors.size() };
        const std::size_t numStrikes{ kMatrix[0].size() };

        UV_REQUIRE(
            wMMatrix[0].size() == numStrikes,
            ErrorCode::InvalidArgument,
            "validateInputs: wMMatrix columns must equal kMatrix columns"
        );

        for (std::size_t i = 1; i < numTenors; ++i)
        {
            UV_REQUIRE(
                tenors[i] > tenors[i - 1],
                ErrorCode::InvalidArgument,
                "validateInputs: tenors must be strictly increasing"
            );
        }

        for (std::size_t i = 0; i < numTenors; ++i)
        {
            UV_REQUIRE(
                kMatrix[i].size() == numStrikes,
                ErrorCode::InvalidArgument,
                "validateInputs: kMatrix is not rectangular"
            );

            UV_REQUIRE(
                wMMatrix[i].size() == numStrikes,
                ErrorCode::InvalidArgument,
                "validateInputs: wMMatrix is not rectangular or does not match kMatrix"
            );

            for (std::size_t j = 1; j < numStrikes; ++j)
            {
                UV_REQUIRE(
                    kMatrix[i][j] > kMatrix[i][j - 1],
                    ErrorCode::InvalidArgument,
                    "validateInputs: kMatrix rows must be strictly increasing"
                );
            }
        }
    }

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

    std::array<double, 4> lowerBounds(const double logKFMin) noexcept
    {
        return
        {
            0.001,                               // b            
            -0.9999,                             // rho
            5.0 * logKFMin,                      // m
            0.01                                 // sigma
        };
    }

    std::array<double, 4> upperBounds(const double logKFMax) noexcept
    {
        return
        {
            2.0,                                // b
            0.9999,                             // rho
            5.0 * logKFMax,                     // m
            10.0                                // sigma
        };
    }

    double calculateWk(const double a,
        const double b,
        const double rho,
        const double m,
        const double sigma,
        const double k) noexcept
    {
        const double x{ k - m };
        return std::fma(b, (rho * x + std::hypot(x, sigma)), a);
    }

    double aParam(const double atmWK,
        const double b,
        const double rho,
        const double m,
        const double sigma) noexcept
    {
        return atmWK - b * (-rho * m + std::sqrt(m * m + sigma * sigma));
    }

    double gk(const GkCache& p) noexcept
    {
        // g(k) = A^2 - B * (w')^2 / 4 + w''/2
        return (p.A * p.A) - 0.25 * p.wkD1Squared * p.B + p.wkD2 * 0.5;  
    }

    std::array<double, 4> gkGrad(
        const double b,
        const double rho,
        const double m,
        const double sigma,
        const double k,
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
        const double wk{ p.wk };
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