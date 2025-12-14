// SPDX-License-Identifier: Apache-2.0
/*
 * File:        SliceData.cpp
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


#include "Math/Functions.hpp"
#include "Core/SliceData.hpp"

#include "Utils/Aux/Errors.hpp"    
#include "Utils/IO/Log.hpp"

#include <cmath>
#include <numbers>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <string> 

namespace uv::core
{
    SliceData::SliceData(Real T,
        const Vector<Real>& mny,
        const Vector<Real>& vol,
        const MarketData& mktData) :
        T_(T),
        mny_(mny),
        vol_(vol),
        mktData_(mktData),
        numStrikes_(mny_.size())
    {
        // Check matching dimensions
        UV_REQUIRE(
            vol_.size() == numStrikes_,
            ErrorCode::InvalidArgument,
            "SliceData::constructor size mismatch: vol_.size() = " + std::to_string(vol_.size()) +
            ", numStrikes_ = " + std::to_string(numStrikes_)
        );

        // Check for positive maturity
        UV_REQUIRE(
            T_ > Real(0.0),
            ErrorCode::InvalidArgument,
            "SliceData::constructor: maturity T_ must be positive (T_ = " + std::to_string(T_) + ")"
        );

        // Check for positive volatility
        auto it = std::find_if(vol_.begin(), vol_.end(),
            [](Real v) { return v <= Real(0.0) || !std::isfinite(v); });
        UV_REQUIRE(
            it == vol_.end(),
            ErrorCode::InvalidArgument,
            "SliceData::constructor: all volatilities must be positive and finite, found invalid value: " +
            std::to_string((it != vol_.end()) ? *it : -Real(1.0))
        );

        // Set variables
        Real r{ mktData_.r };
        Real q{ mktData_.q };
        Real S{ mktData_.S };

        // Calculate forward price
        F_ = S * std::exp((r - q) * T_);

        // Set the strikes vector
        K_.resize(numStrikes_);
        std::transform(mny_.begin(), mny_.end(), K_.begin(),
            [S](Real mny)
            {
                return mny * S;
            });

        // Convert plain strikes K into log-moneyness log(K/F)
        logKF_.resize(numStrikes_);
        std::transform(K_.begin(), K_.end(), logKF_.begin(),
            [F = F_](Real K)
            {
                return std::log(K / F);
            });

        // Convert implied volatility into total variance
        wT_.resize(numStrikes_);
        std::transform(vol_.begin(), vol_.end(), wT_.begin(),
            [T = T_](Real vol)
            {
                return vol * vol * T;
            });

        // Calculate Call Black-Scholes price
        callBS_.resize(numStrikes_);
        for (std::size_t i = 0; i < numStrikes_; ++i)
        {
            callBS_[i] = math::blackScholes
            (
                T_,
                r,
                q,
                vol_[i],
                S,
                K_[i],
                /*isCall=*/true
            );
        }
    }

    Real SliceData::T() const noexcept
    {
        return T_;
    }

    Real SliceData::F() const noexcept
    {
        return F_;
    }

    std::size_t SliceData::numStrikes() const noexcept
    {
        return numStrikes_;
    }

    Real SliceData::r() const noexcept
    {
        return mktData_.r;
    }

    const Vector<Real>& SliceData::mny() const noexcept
    {
        return mny_;
    }

    const Vector<Real>& SliceData::logKF() const noexcept
    {
        return logKF_;
    }

    const Vector<Real>& SliceData::vol() const noexcept
    {
        return vol_;
    }

    const Vector<Real>& SliceData::wT() const noexcept
    {
        return wT_;
    }

    const Vector<Real>& SliceData::K() const noexcept
    {
        return K_;
    }

    const Vector<Real>& SliceData::callBS() const noexcept
    {
        return callBS_;
    }

    void SliceData::setWt(const Vector<Real>& wT)
    {
        // Check matching size
        UV_REQUIRE(wT.size() == numStrikes_, ErrorCode::InvalidArgument,
            "SliceData::setWT size mismatch: got " + std::to_string(wT.size()) +
            ", expected " + std::to_string(numStrikes_));

        // Validate all total variances are non-negative before mutating state
        for (std::size_t i = 0; i < numStrikes_; ++i)
        {
            UV_REQUIRE(wT[i] >= Real(0.0), ErrorCode::InvalidArgument,
                "SliceData::setWT: negative total variance at index " + std::to_string(i) +
                " (value = " + std::to_string(wT[i]) + ")");
        }

        // Create a copy
        wT_ = wT;

        // Recalculate volatility per strike: sigma = sqrt(w/T)
        for (size_t i = 0; i < numStrikes_; ++i)
        {
            vol_[i] = std::sqrt(wT_[i] / T_);
        }

        // Recalculate Call Black-Scholes prices
        Real r{ mktData_.r };
        Real q{ mktData_.q };
        Real S{ mktData_.S };

        for (std::size_t i = 0; i < numStrikes_; ++i)
        {
            callBS_[i] = math::blackScholes
            (
                T_,
                r,
                q,
                vol_[i],
                S,
                K_[i],
                /*isCall=*/true
            );
        }
    }

    void SliceData::setCallBS(const Vector<Real>& callBS)
    {
        // Check matching size
        UV_REQUIRE(callBS.size() == numStrikes_, ErrorCode::InvalidArgument,
            "SliceData::setCallBS size mismatch: got " + std::to_string(callBS.size()) +
            ", expected " + std::to_string(numStrikes_));

        // All option prices must be positive
        for (std::size_t i = 0; i < callBS.size(); ++i)
        {
            UV_REQUIRE(callBS[i] > Real(0.0), ErrorCode::InvalidArgument,
                "SliceData::setCallBS: negative or zero call price at index " + std::to_string(i) +
                " (value = " + std::to_string(callBS[i]) + ")");
        }

        // Create a copy
        callBS_ = callBS;

        // Recalculate volatility and total variance
        for (size_t i = 0; i < numStrikes_; ++i)
        {   
            Real vol
            {
                math::impliedVolBS
                (
                    callBS_[i],
                    T_,
                    mktData_.r,
                    mktData_.q,
                    mktData_.S,
                    K_[i],
                    /*isCall=*/true
                    // ftolAbs & maxEval use defaults from API
                )
            };

            vol_[i] = vol;
            wT_[i] = vol * vol * T_;
        }
    }
}
