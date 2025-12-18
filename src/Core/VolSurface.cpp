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
#include "Core/Matrix/Functions.hpp"
#include "Core/MarketData.hpp"
#include "Math/Functions.hpp"
#include "Utils/IO/Functions.hpp"
#include "Utils/Types.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <span>

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
        rates_(numTenors_, mktData.r),          // Constant risk-free rate for now
        dividends_(numTenors_, mktData.q),      // Constant dividends for now
        strikes_(numStrikes_),
        forwards_(numTenors_),
        callPrices_(numTenors_, numStrikes_),
        logKFMatrix_(numTenors_, numStrikes_),
        totVarMatrix_(numTenors_, numStrikes_)
    {   

        // ---------- Initialize member variables ---------- 

        strikes_ = core::multiply(mny_, S_);
        setForwards_();
        setCallPrices_();
        setLogKFMatrix_();
        setTotVar_(volMatrix_);
    }

    void VolSurface::setForwards_() noexcept
    {
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            forwards_[i] = S_ * std::exp((rates_[i] - dividends_[i]) * tenors_[i]);
        }
    }

    void VolSurface::setCallPrices_()
    {
        callPrices_ = applyIndexed
            (
                volMatrix_,
                [&](std::size_t i, std::size_t j, Real sigma)
                {
                    return math::blackScholes(
                        tenors_[i],
                        rates_[i],
                        dividends_[i],
                        sigma,
                        S_,
                        strikes_[j]
                    );
                }
            );
    }

    void VolSurface::setLogKFMatrix_()
    {
        logKFMatrix_ = applyIndexed<Real>
            (
                numTenors_,
                numStrikes_,
                [&](std::size_t i, std::size_t j)
                {
                    return std::log(strikes_[j]) - std::log(forwards_[i]);
                }
            );
    }

    void VolSurface::setTotVar_(const Matrix<Real>& volMatrix)
    {
        totVarMatrix_ = volMatrix;
        squareInplace(totVarMatrix_);
        hadamardInplace(totVarMatrix_, tenors_);
    }

    void VolSurface::setVolFromVar_(const Matrix<Real>& totVarMatrix)
    {
        volMatrix_ = totVarMatrix;
        hadamardInplace(volMatrix_, reciprocal(tenors_));
        sqrtInplace(volMatrix_);
    }

    void VolSurface::setVolFromPrices_(const Matrix<Real>& callPrices)
    {
        volMatrix_ = applyIndexed
        (
            callPrices,
            [&](std::size_t i, std::size_t j, Real callPrice)
            {
                return math::impliedVolBS(
                    callPrice,
                    tenors_[i],
                    rates_[i],
                    dividends_[i],
                    S_,
                    strikes_[j]
                );
            }
        );
    }

    void VolSurface::setTotVar(const Matrix<Real>& totVarMatrix) 
    {
        UV_REQUIRE
        (
            totVarMatrix.rows() == numTenors_ &&
            totVarMatrix.cols() == numStrikes_,
            ErrorCode::InvalidArgument,
            "setTotVar: dimension mismatch"
        );

        totVarMatrix_ = totVarMatrix;
        setVolFromVar_(totVarMatrix);
        setCallPrices_();
    }

    void VolSurface::setCallPrices(const Matrix<Real>& callPrices)
    {
        UV_REQUIRE
        (
            callPrices.rows() == numTenors_ &&
            callPrices.cols() == numStrikes_,
            ErrorCode::InvalidArgument,
            "setCallPrices: dimension mismatch"
        );

        callPrices_ = callPrices;
        setVolFromPrices_(callPrices_);
        setTotVar_(volMatrix_);
    }

    const Vector<Real>& VolSurface::tenors() const noexcept
    {
        return tenors_;
    }

    std::size_t VolSurface::numTenors() const noexcept
    {
        return numTenors_;
    }

    const Vector<Real>& VolSurface::mny() const noexcept
    {
        return mny_;
    }

    std::size_t VolSurface::numStrikes() const noexcept
    {
        return numStrikes_;
    }

    const Matrix<Real>& VolSurface::volMatrix() const noexcept
    {
        return volMatrix_;
    }

    Real VolSurface::S() const noexcept
    {
        return S_;
    }

    const Vector<Real>& VolSurface::rates() const noexcept
    {
        return rates_;
    }

    const Vector<Real>& VolSurface::dividends() const noexcept
    {
        return dividends_;
    }

    const Vector<Real>& VolSurface::strikes() const noexcept
    {
        return strikes_;
    }

    const Vector<Real>& VolSurface::forwards() const noexcept
    {
        return forwards_;
    }

    const Matrix<Real>& VolSurface::callPrices() const noexcept
    {
        return callPrices_;
    }

    const Matrix<Real>& VolSurface::logKFMatrix() const noexcept
    {
        return logKFMatrix_;
    }

    const Matrix<Real>& VolSurface::totVarMatrix() const noexcept
    {
        return totVarMatrix_;
    }

    void VolSurface::printVol(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const Real> header
        {
            mnyFlag ? std::span<const Real>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\log(F/K)" };

        utils::printMatrix(
            /*title=*/title,
            /*header=*/header,
            /*rowLabels=*/tenors_,
            /*M=*/volMatrix_,
            /*headerPrec=*/headerPrec,
            /*rowLabelPrec=*/2,
            /*valuePrec=*/valuePrec
        );
    }

    void VolSurface::printTotVar(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const Real> header
        { 
            mnyFlag ? std::span<const Real>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\log(F/K)" };

        utils::printMatrix(
            /*title=*/title,
            /*header=*/header,
            /*rowLabels=*/tenors_,
            /*M=*/totVarMatrix_,
            /*headerPrec=*/headerPrec,
            /*rowLabelPrec=*/2,
            /*valuePrec=*/valuePrec
        );
    }

    void VolSurface::printBSCall(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const Real> header
        {
            mnyFlag ? std::span<const Real>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\log(F/K)" };

        utils::printMatrix(
            /*title=*/title,
            /*header=*/header,
            /*rowLabels=*/tenors_,
            /*M=*/callPrices_,
            /*headerPrec=*/headerPrec,
            /*rowLabelPrec=*/2,
            /*valuePrec=*/valuePrec
        );
    }
}
