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

        Vector<Real> tenors_;                       
        std::size_t numTenors_;                     

        Vector<Real> mny_;                         
        std::size_t numStrikes_;

        Matrix<Real> volMatrix_;

        Real S_;
        Vector<Real> rates_;
        Vector<Real> dividends_;

        Vector<Real> strikes_;
        Vector<Real> forwards_;

        Matrix<Real> callPrices_;
        Matrix<Real> logKFMatrix_; 
        Matrix<Real> totVarMatrix_;


    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------

        VolSurface() = delete;
        explicit VolSurface(Vector<Real> tenors,
            Vector<Real> mny,
            Matrix<Real> volMatrix,
            const MarketData& mktData);

        // Print volatility surface on the console
        void printVol() const noexcept;

        // Print total variance surface on the console
        void printTotVar() const noexcept;

        // Print Black-Scholes Call price surface
        void printBSCall() const noexcept;


        // Getters
        const Matrix<Real>& callPrices() const noexcept;
        const Matrix<Real>& logKFMatrix() const noexcept;
        const Matrix<Real>& totVarMatrix() const noexcept;
        const Vector<Real>& forwards() const noexcept;
        const Matrix<Real>& volMatrix() const noexcept;
        const Vector<Real>& tenors() const noexcept;
        const Vector<Real>& strikes() const noexcept;
        std::size_t numTenors() const noexcept;
        std::size_t numStrikes() const noexcept;
        const Vector<Real>& rates() const noexcept;
        const Vector<Real>& dividends() const noexcept;

        // Setters
        void setTotVar(const Matrix<Real>& totalVarMatrix);
        void setCallPrices(const Matrix<Real>& callPrices);
    };
}
