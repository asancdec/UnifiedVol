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


#include "Core/Functions.hpp"
#include "Core/Matrix/Functions.hpp"
#include "Core/MarketData.hpp"
#include "Math/Functions.hpp"
#include "Utils/IO/Functions.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <span>

namespace uv::core
{
    template <std::floating_point T>
    VolSurface<T>::VolSurface(Vector<T> tenors,
        Vector<T> mny,
        Matrix<T> volMatrix,
        const MarketData<T>& mktData)
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
        discountFactors_(numTenors_),
        callPrices_(numTenors_, numStrikes_),
        logKFMatrix_(numTenors_, numStrikes_),
        totVarMatrix_(numTenors_, numStrikes_)
    {

        // ---------- Initialize member variables ---------- 

        strikes_ = core::multiply<T>(mny_, S_);
        setDiscountFactors_();
        setForwards_();
        setCallPrices_();
        setLogKFMatrix_();
        setTotVar_(volMatrix_);
    }

    template <std::floating_point T>
    void VolSurface<T>::setDiscountFactors_() noexcept
    {
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            discountFactors_[i] = std::exp(-rates_[i] * tenors_[i]);
        }
    }

    template <std::floating_point T>
    void VolSurface<T>::setForwards_() noexcept
    {
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            forwards_[i] = S_ * std::exp((rates_[i] - dividends_[i]) * tenors_[i]);
        }
    }

    template <std::floating_point T>
    void VolSurface<T>::setCallPrices_()
    {
        callPrices_ = transformIndexed<T>
            (
                volMatrix_,
                [&](std::size_t i, std::size_t j, T sigma)
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

    template <std::floating_point T>
    void VolSurface<T>::setLogKFMatrix_()
    {
        logKFMatrix_ = generateIndexed<T>
            (
                numTenors_,
                numStrikes_,
                [&](std::size_t i, std::size_t j)
                {
                    return std::log(strikes_[j]) - std::log(forwards_[i]);
                }
            );
    }

    template <std::floating_point T>
    void VolSurface<T>::setTotVar_(const Matrix<T>& volMatrix)
    {
        totVarMatrix_ = volMatrix;
        squareInplace(totVarMatrix_);
        hadamardInplace(totVarMatrix_, tenors_);
    }

    template <std::floating_point T>
    void VolSurface<T>::setVolFromVar_(const Matrix<T>& totVarMatrix)
    {
        volMatrix_ = totVarMatrix;
        hadamardInplace(volMatrix_, reciprocal<T>(tenors_));
        sqrtInplace(volMatrix_);
    }

    template <std::floating_point T>
    void VolSurface<T>::setVolFromPrices_(const Matrix<T>& callPrices)
    {
        volMatrix_ = transformIndexed<T>
            (
                callPrices,
                [&](std::size_t i, std::size_t j, T callPrice)
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

    template <std::floating_point T>
    void VolSurface<T>::setTotVar(const Matrix<T>& totVarMatrix)
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

    template <std::floating_point T>
    void VolSurface<T>::setCallPrices(const Matrix<T>& callPrices)
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

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::tenors() const noexcept
    {
        return tenors_;
    }

    template <std::floating_point T>
    std::size_t VolSurface<T>::numTenors() const noexcept
    {
        return numTenors_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::mny() const noexcept
    {
        return mny_;
    }

    template <std::floating_point T>
    std::size_t VolSurface<T>::numStrikes() const noexcept
    {
        return numStrikes_;
    }

    template <std::floating_point T>
    const Matrix<T>& VolSurface<T>::volMatrix() const noexcept
    {
        return volMatrix_;
    }

    template <std::floating_point T>
    T VolSurface<T>::S() const noexcept
    {
        return S_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::rates() const noexcept
    {
        return rates_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::dividends() const noexcept
    {
        return dividends_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::discountFactors() const noexcept
    {
        return discountFactors_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::strikes() const noexcept
    {
        return strikes_;
    }

    template <std::floating_point T>
    const Vector<T>& VolSurface<T>::forwards() const noexcept
    {
        return forwards_;
    }

    template <std::floating_point T>
    const Matrix<T>& VolSurface<T>::callPrices() const noexcept
    {
        return callPrices_;
    }

    template <std::floating_point T>
    const Matrix<T>& VolSurface<T>::logKFMatrix() const noexcept
    {
        return logKFMatrix_;
    }

    template <std::floating_point T>
    const Matrix<T>& VolSurface<T>::totVarMatrix() const noexcept
    {
        return totVarMatrix_;
    }

    template <std::floating_point T>
    void VolSurface<T>::printVol(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const T> header
        {
            mnyFlag ? std::span<const T>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\lnKF" };

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

    template <std::floating_point T>
    void VolSurface<T>::printTotVar(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const T> header
        {
            mnyFlag ? std::span<const T>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\lnKF" };

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

    template <std::floating_point T>
    void VolSurface<T>::printBSCall(const unsigned int valuePrec,
        const bool mnyFlag) const noexcept
    {
        const std::span<const T> header
        {
            mnyFlag ? std::span<const T>(mny_) : logKFMatrix_[0]
        };
        const int headerPrec{ mnyFlag ? 2 : 4 };
        const char* title{ mnyFlag ? "T\\%S" : "T\\lnKF" };

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
