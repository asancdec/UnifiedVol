/**
* SliceData.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/MathFunctions/MathFunctions.hpp"
#include "Core/SliceData.hpp"
#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"
#include <cmath>
#include <numbers>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <string> 

namespace uv
{
    SliceData::SliceData(double T,
        const ::std::vector<double>& mny,
        const ::std::vector<double>& vol,
        const MarketData& mktData) :
        T_(T),
        mny_(mny),
        vol_(vol),
        mktData_(mktData)
    {
        // Check matching dimensions
        UV_REQUIRE(
            vol_.size() == mny_.size(),
            ErrorCode::InvalidArgument,
            "SliceData::constructor size mismatch: vol_.size() = " + ::std::to_string(vol_.size()) +
            ", mny_.size() = " + ::std::to_string(mny_.size())
        );

        // Check for positive maturity
        UV_REQUIRE(
            T_ > 0.0,
            ErrorCode::InvalidArgument,
            "SliceData::constructor: maturity T_ must be positive (T_ = " + ::std::to_string(T_) + ")"
        );

        // Check for positive volatility
        auto it = ::std::find_if(vol_.begin(), vol_.end(),
            [](double v) { return v <= 0.0 || !::std::isfinite(v); });
        UV_REQUIRE(
            it == vol_.end(),
            ErrorCode::InvalidArgument,
            "SliceData::constructor: all volatilities must be positive and finite, found invalid value: " +
            ::std::to_string((it != vol_.end()) ? *it : -1.0)
        );

        // Set variables
        double r{ mktData_.r };
        double q{ mktData_.q };
        double S{ mktData_.S };

        // Calculate forward price
        F_ = S * ::std::exp((r - q) * T_);

        // Set the strikes vector
        K_.resize(mny_.size());
        ::std::transform(mny_.begin(), mny_.end(), K_.begin(),
            [S](double mny)
            {
                return mny * S;
            });

        // Convert plain strikes K into log-moneyness log(K/F)
        logFM_.resize(K_.size());
        ::std::transform(K_.begin(), K_.end(), logFM_.begin(),
            [F = F_](double K)
            {
                return ::std::log(K / F);
            });

        // Convert implied volatility into total variance
        wT_.resize(vol_.size());
        ::std::transform(vol_.begin(), vol_.end(), wT_.begin(),
            [T = T_](double vol)
            {
                return vol * vol * T;
            });

        // Calculate Call Black-Scholes price
        callBS_.resize(vol_.size());
        for (::std::size_t i = 0; i < vol_.size(); ++i)
        {
            callBS_[i] = blackScholes
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

    double SliceData::minWT() const noexcept
    {
        return *::std::min_element(wT_.begin(), wT_.end());
    }

    double SliceData::maxWT() const noexcept
    {
        return *::std::max_element(wT_.begin(), wT_.end());
    }

    double SliceData::minLogFM() const noexcept
    {
        return *::std::min_element(logFM_.begin(), logFM_.end());
    }

    double SliceData::maxLogFM() const noexcept
    {
        return *::std::max_element(logFM_.begin(), logFM_.end());
    }

    double SliceData::atmWT() const noexcept
    {
        const ::std::size_t n{ (::std::min)(logFM_.size(), wT_.size()) };
        const double* first{ logFM_.data() };
        const double* last{ first + n };
        const double* it = ::std::min_element(first, last, [](double a, double b)
            {
                return ::std::abs(a) < ::std::abs(b);
            });
        return wT_[static_cast<::std::size_t>(it - first)];
    }

    double SliceData::T() const noexcept
    {
        return T_;
    }

    double SliceData::F() const noexcept
    {
        return F_;
    }

    double SliceData::r() const noexcept
    {
        return mktData_.r;
    }

    const ::std::vector<double>& SliceData::mny() const noexcept
    {
        return mny_;
    }

    const ::std::vector<double>& SliceData::logFM() const noexcept
    {
        return logFM_;
    }

    const ::std::vector<double>& SliceData::vol() const noexcept
    {
        return vol_;
    }

    const ::std::vector<double>& SliceData::wT() const noexcept
    {
        return wT_;
    }

    const ::std::vector<double>& SliceData::K() const noexcept
    {
        return K_;
    }

    const ::std::vector<double>& SliceData::callBS() const noexcept
    {
        return callBS_;
    }

    void SliceData::setWT(const ::std::vector<double>& wT)
    {
        // Check matching size
        UV_REQUIRE(wT.size() == wT_.size(), ErrorCode::InvalidArgument,
            "SliceData::setWT size mismatch: got " + ::std::to_string(wT.size()) +
            ", expected " + ::std::to_string(wT_.size()));

        // Validate all total variances are non-negative before mutating state
        for (std::size_t i = 0; i < wT.size(); ++i)
        {
            UV_REQUIRE(wT[i] >= 0.0, ErrorCode::InvalidArgument,
                "SliceData::setWT: negative total variance at index " + ::std::to_string(i) +
                " (value = " + ::std::to_string(wT[i]) + ")");
        }

        // Create a copy
        wT_ = wT;

        // Recalculate volatility per strike: sigma = sqrt(w/T)
        const double invT{ 1.0 / T_ };
        for (size_t i = 0; i < wT_.size(); ++i)
        {
            vol_[i] = ::std::sqrt(wT_[i] * invT);
        }

        // Recalculate Call Black-Scholes prices
        double r{ mktData_.r };
        double q{ mktData_.q };
        double S{ mktData_.S };

        for (::std::size_t i = 0; i < vol_.size(); ++i)
        {
            callBS_[i] = blackScholes
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

    void SliceData::setCallBS(const ::std::vector<double>& callBS)
    {
        // Check matching size
        UV_REQUIRE(callBS.size() == callBS_.size(), ErrorCode::InvalidArgument,
            "SliceData::setCallBS size mismatch: got " + ::std::to_string(callBS.size()) +
            ", expected " + ::std::to_string(callBS_.size()));

        // All option prices must be positive
        for (std::size_t i = 0; i < callBS.size(); ++i)
        {
            UV_REQUIRE(callBS[i] > 0.0, ErrorCode::InvalidArgument,
                "SliceData::setCallBS: negative or zero call price at index " + ::std::to_string(i) +
                " (value = " + ::std::to_string(callBS[i]) + ")");
        }

        // Create a copy
        callBS_ = callBS;

        // Recalculate volatility and total variance
        for (size_t i = 0; i < callBS_.size(); ++i)
        {   
            double vol
            {
                uv::impliedVolBS
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
