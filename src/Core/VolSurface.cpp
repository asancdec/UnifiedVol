// SPDX-License-Identifier: Apache-2.0
/*
 * File:        VolSurface.cpp
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


#include "Core/VolSurface.hpp"
#include "Core/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"
#include "Math/Functions.hpp"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>

namespace uv::core
{
    VolSurface::VolSurface(Vector<Real> tenors,
        Vector<Real> mny,
        Matrix<Real> volMatrix,
        const MarketData& mktData)
        : tenors_(std::move(tenors)),
        numTenors_(tenors_.size()),
        mny_(std::move(mny)),
        numStrikes_(mny_.size()),
        volMatrix_(std::move(volMatrix)),
        S_(mktData.S),
        rates_(numTenors_, mktData.r),         // Constant risk-free rate for now
        dividends_(numTenors_, mktData.q),     // Constant dividends for now
        strikes_(uv::core::multiply(mny_, S_)),
        forwards_(numTenors_),
        callPrices_(uv::math::blackScholes(tenors_, rates_, dividends_, volMatrix_, S_, strikes_)),
        logKFMatrix_(numTenors_, Vector<Real>(numStrikes_)),
        totVarMatrix_(numTenors_, Vector<Real>(numStrikes_))
    {   

        // ---------- Initialize member variables ---------- 

        // Forwards
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            forwards_[i] = S_ * std::exp((rates_[i] - dividends_[i]) * tenors_[i]);
        }

        // LogKFMatrix
        for (std::size_t i = 0; i < numTenors_; ++i)
        {   
            const Real logF{ std::log(forwards_[i]) };
            Vector<Real>& row{ logKFMatrix_[i] };

            for (std::size_t j = 0; j < numStrikes_; ++j)
            {
                row[j] = std::log(strikes_[j]) - logF;
            }
        }

        // TotVarMatrix
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            totVarMatrix_[i] = uv::core::multiply(uv::core::hadamard(volMatrix_[i], volMatrix_[i]), tenors_[i]);
        }
    }

    void VolSurface::printVol() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\%S\t";

        // Header row (moneyness)
        for (const auto& m : mny_)
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : volMatrix_[i])
                oss << std::fixed << std::setprecision(5) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printTotVar() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : mny_)
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : totVarMatrix_[i])

                oss << std::fixed << std::setprecision(5) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printBSCall() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : mny_)
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : callPrices_[i])
                oss << std::fixed << std::setprecision(5) << v << '\t';

            oss << '\n';
        }
        UV_INFO(oss.str());
    }

    const Matrix<Real>& VolSurface::totVarMatrix() const noexcept
    {
      
        return totVarMatrix_;
    }

    const Matrix<Real>& VolSurface::volMatrix() const noexcept
    {      
        return volMatrix_;
    }

    const Matrix<Real>& VolSurface::logKFMatrix() const noexcept
    {
        return logKFMatrix_;
    }

    const Vector<Real>& VolSurface::forwards() const noexcept
    {
        return forwards_;
    }

    const Matrix<Real>& VolSurface::callPrices() const noexcept
    {
        return callPrices_;
    }

    const Vector<Real>& VolSurface::tenors() const noexcept
    {
        return tenors_;
    }

    const Vector<Real>& VolSurface::strikes() const noexcept
    {
        return strikes_;
    }

    std::size_t VolSurface::numTenors() const noexcept
    {
        return numTenors_;
    }

    std::size_t VolSurface::numStrikes() const noexcept
    {
        return numStrikes_;
    }

    const Vector<Real>& VolSurface::rates() const noexcept
    {
        return rates_;
    }

    const Vector<Real>& VolSurface::dividends() const noexcept
    {
        return dividends_;
    }

    void VolSurface::setTotVar(const Matrix<Real>& totVarMatrix)
    {
        // ---------- Set ----------

        totVarMatrix_ = totVarMatrix;

        for (std::size_t i = 0; i < totVarMatrix_.size(); ++i)
        {
            const Real invT = Real(1) / tenors_[i];

            for (std::size_t j = 0; j < totVarMatrix_[i].size(); ++j)
            {
                volMatrix_[i][j] =
                    std::sqrt(totVarMatrix_[i][j] * invT);
            }
        }

        callPrices_ = uv::math::blackScholes(tenors_, rates_, dividends_, volMatrix_, S_, strikes_);

    }

    void VolSurface::setCallPrices(const Matrix<Real>& callPrices)
    {
        // ---------- Set ----------

        callPrices_ = callPrices;


        for (std::size_t i = 0; i < callPrices_.size(); ++i)
        {
            const Real T = tenors_[i];
            const Real r = rates_[i];
            const Real q = dividends_[i];

            for (std::size_t j = 0; j < callPrices_[i].size(); ++j)
            {
                volMatrix_[i][j] = uv::math::impliedVolBS(
                    callPrices_[i][j],
                    T,
                    r,
                    q,
                    S_,
                    strikes_[j],
                    /*isCall=*/true
                );
            }
        }


    }
}
