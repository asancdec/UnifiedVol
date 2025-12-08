// SPDX-License-Identifier: Apache-2.0
/*
 * File:        SliceData.hpp
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

#include "Core/MarketData.hpp"
#include "Utils/Types.hpp"

namespace uv::core
{
    class SliceData
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        Real T_;                    // Slice maturity
        Real F_;                    // Forward price
        Vector<Real> mny_;          // Plain moneyness (K/S)
        Vector<Real> logFM_;        // Log-forward moneyness log(K/F)
        Vector<Real> vol_;          // Volatilities 
        Vector<Real> wT_;           // Total variance (vol²T) 
        Vector<Real> K_;            // Strikes vector
        Vector<Real> callBS_;       // Value of European call options
        MarketData mktData_;        // Market data struct
        std::size_t numStrikes_;    // Number of strikes

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------   
        SliceData() = delete;
        explicit SliceData(Real T,
            const Vector<Real>& mny,
            const Vector<Real>& vol,
            const MarketData& mktData);

        //--------------------------------------------------------------------------
        // Math functions
        //--------------------------------------------------------------------------

        // Determine minimum total variance
        Real minWT() const noexcept;

        // Determine maximum total variance
        Real maxWT() const noexcept;

        // Determine minimum log-forward moneyness
        Real minLogFM() const noexcept;

        // Determine maximum log-forward moneyness
        Real maxLogFM() const noexcept;

        // Return ATM total variance
        Real atmWT() const noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        Real T() const noexcept;
        Real F() const noexcept;
        std::size_t numStrikes() const noexcept;
        Real r() const noexcept;
        const Vector<Real>& mny() const noexcept;
        const Vector<Real>& logFM() const noexcept;
        const Vector<Real>& vol() const noexcept;
        const Vector<Real>& wT() const noexcept;
        const Vector<Real>& K() const noexcept;
        const Vector<Real>& callBS() const noexcept;


        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------
        void setWT(const Vector<Real>& wT);
        void setCallBS(const Vector<Real>& callBS);
    };
}
