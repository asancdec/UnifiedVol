// SPDX-License-Identifier: Apache-2.0
/*
 * File:        VolSurface.hpp
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


#pragma once

#include "Core/SliceData.hpp"
#include "Core/MarketData.hpp"
#include "Utils/Types.hpp"

#include <cstddef>

namespace uv::core
{
    class VolSurface
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        std::vector<SliceData> slices_;         // 2D volatility grid: vols[maturity][tenors slice]
        Vector<Real> tenors_;                   // Vector with different volatility surface tenors
        Vector<Real> strikes_;                  // Strikes

        std::size_t numTenors_;                 // Number of tenors in years
        std::size_t numStrikes_;                // Number of strikes

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------

        VolSurface() = delete;
        explicit VolSurface(const Vector<Real>& mny,
            const Matrix<Real>& vols,
            const Vector<Real>& tenors,
            const MarketData& mktData);

        //--------------------------------------------------------------------------
        // Utilities
        //--------------------------------------------------------------------------

        // Print volatility surface on the console
        void printVol() const noexcept;

        // Print total variance surface on the console
        void printTotVar() const noexcept;

        // Print Black-Scholes Call price surface
        void printBSCall() const noexcept;

        // Return total variance matrix
        Matrix<Real> totVarMatrix() const noexcept;

        // Return total variance matrix
        Matrix<Real> volMatrix() const noexcept;

        // Return variance matrix
        Matrix<Real> varMatrix() const noexcept;

        // Return the logKF matrix
        Matrix<Real> logKFMatrix() const noexcept;

        Vector<Real> forwards() const noexcept;
        Matrix<Real> calls() const noexcept;

        // Getters
        std::vector<SliceData>& slices() noexcept;
        const Vector<Real>& tenors() const noexcept;
        const Vector<Real>& strikes() const noexcept;
        std::size_t numTenors() const noexcept;
        std::size_t numStrikes() const noexcept;
        Vector<Real> rates() const noexcept;

        // Setters
        void setWt(const Matrix<Real>& wT);
        void setCallBS(const Matrix<Real>& calls);
    };
}
