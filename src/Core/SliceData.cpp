/**
* SliceData.cpp
* Author: Alvaro Sanchez de Carlos
*/

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

using uv::ErrorCode;

SliceData::SliceData(double T,
    const std::vector<double>& mny,
    const std::vector<double>& vol,
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
        "SliceData::constructor size mismatch: vol_.size() = " + std::to_string(vol_.size()) +
        ", mny_.size() = " + std::to_string(mny_.size())
    );

    // Check for positive maturity
    UV_REQUIRE(
        T_ > 0.0,
        ErrorCode::InvalidArgument,
        "SliceData::constructor: maturity T_ must be positive (T_ = " + std::to_string(T_) + ")"
    );

    // Check for positive volatility
    auto it = std::find_if(vol_.begin(), vol_.end(),
        [](double v) { return v <= 0.0 || !std::isfinite(v); });
    UV_REQUIRE(
        it == vol_.end(),
        ErrorCode::InvalidArgument,
        "SliceData::constructor: all volatilities must be positive and finite, found invalid value: " +
        std::to_string((it != vol_.end()) ? *it : -1.0)
    );

    // Set variables
    double r{ mktData_.r };
    double q{ mktData_.q };
    double S{ mktData_.S };

    // Calculate forward price
    F_ = S * std::exp((r - q) * T_);

    // Set the strikes vector
    K_.resize(mny_.size());
    std::transform(mny_.begin(), mny_.end(), K_.begin(), 
        [S](double mny)
        {
            return mny * S;
        });

    // Convert plain strikes K into log-moneyness log(K/F)
    logFM_.resize(K_.size());
    std::transform(K_.begin(), K_.end(), logFM_.begin(), 
        [F = F_](double K)
        {       
            return std::log(K / F); 
        });

    // Convert implied volatility into total variance
    wT_.resize(vol_.size());
    std::transform(vol_.begin(), vol_.end(),wT_.begin(),
        [T = T_](double vol)
        {
            return vol * vol * T;
        });

    // Calculate Call Black-Scholes price
    callBS_.resize(vol_.size());
    for (std::size_t i = 0; i < vol_.size(); ++i)
    {
        callBS_[i] = blackScholes
        (
            T_,
            r,
            q,
            vol_[i],
            S,
            K_[i],
            true
        );
    }
}

double SliceData::blackScholes(double T,
    double r,
    double q,
    double vol,
    double S,
    double K,
    bool isCall) noexcept
{
    double volSqrtT{ std::sqrt(T) * vol };
    double d1{std::fma(T, (r - q + vol * vol * 0.5), std::log(S/K)) / volSqrtT };
    double d2{ std::fma(-1.0, volSqrtT, d1) };

    if (isCall)
    {
        return S * std::exp(-q * T) * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
    }
    else
    {
        return K * std::exp(-r * T) * normalCDF(-d2) - S * std::exp(-q * T) * normalCDF(-d1);
    }
}

double SliceData::normalCDF(double x) noexcept
{
    return std::erfc( -x / std::sqrt(2.0)) * 0.5;
}

double SliceData::minWT() const noexcept
{
    return *std::min_element(wT_.begin(), wT_.end());
}

double SliceData::maxWT() const noexcept
{
    return *std::max_element(wT_.begin(), wT_.end());
}

double SliceData::minLogFM() const noexcept
{
    return *std::min_element(logFM_.begin(), logFM_.end());
}

double SliceData::maxLogFM() const noexcept
{
    return *std::max_element(logFM_.begin(), logFM_.end());
}

double SliceData::atmWT() const noexcept
{   
    const std::size_t n{ (std::min)(logFM_.size(), wT_.size()) };
    const double* first{ logFM_.data() };
    const double* last{ first + n };
    const double* it = std::min_element(first, last, [](double a, double b)
        { 
            return std::abs(a) < std::abs(b); 
        });
    return wT_[static_cast<std::size_t>(it - first)];
}

double SliceData::T() const noexcept
{
    return T_;
}

double SliceData::F() const noexcept
{
    return F_;
}

const std::vector<double>& SliceData::mny() const noexcept
{
    return mny_;
}

const std::vector<double>& SliceData::logFM() const noexcept
{
    return logFM_;
}

const std::vector<double>& SliceData::vol() const noexcept
{
    return vol_;
}

const std::vector<double>& SliceData::wT() const noexcept
{
    return wT_;
}

const std::vector<double>& SliceData::K() const noexcept
{
    return K_;
}

const std::vector<double>& SliceData::callBS() const noexcept
{
    return callBS_;
}

void SliceData::setWT(const std::vector<double>& wT)
{   
    // Check matching size
    UV_REQUIRE(wT.size() == wT_.size(), ErrorCode::InvalidArgument,
        "SliceData::setWT size mismatch: got " + std::to_string(wT.size()) +
        ", expected " + std::to_string(wT_.size()));

    // Create a copy
    wT_ = wT;  

    // Recalculate volatility per strike: sigma = sqrt(w/T)
    const double invT{ 1.0 / T_ };
    for (size_t i = 0; i < wT_.size(); ++i) 
    {
        const double wk = wT_[i];
        UV_REQUIRE(wk >= 0.0, ErrorCode::InvalidArgument,
            "Negative total variance at idx " + std::to_string(i));
        vol_[i] = std::sqrt(wk * invT);
    }

    // Recalculate Call Black-Scholes prices
    double r{ mktData_.r };
    double q{ mktData_.q };
    double S{ mktData_.S };

    for (std::size_t i = 0; i < vol_.size(); ++i)
    {
        callBS_[i] = blackScholes
        (
            T_,
            r,
            q,
            vol_[i],
            S,
            K_[i],
            true
        );
    }
}