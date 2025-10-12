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
    const std::vector<double>& logFM,
    const std::vector<double>& vol, 
    const std::vector<double>& wT,
    const MarketData& marketData)
    : T_(T), mny_(mny), logFM_(logFM), vol_(vol), wT_(wT), marketData_(marketData)
{   
    // Set variables
    double r{ marketData_.r };
    double q{ marketData_.q };
    double S{ marketData_.S };

    // Set the strikes vector
    K_.reserve(mny_.size());
    std::transform(mny_.begin(), mny_.end(),
        std::back_inserter(K_), [S](double mny)
        {
            return mny * S;
        });

    // Check matching size
    UV_REQUIRE(
        vol_.size() == K_.size(),
        ErrorCode::InvalidArgument,
        "SliceData::constructor size mismatch: vol_.size() = " + std::to_string(vol_.size()) +
        ", K_.size() = " + std::to_string(K_.size())
    );

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

SliceData SliceData::fromMarketData(const std::vector<double>& mny,
    const std::vector<double>& vol,
    const MarketData& mkt, double T)
{   
    // Compute forward factor F/S0
    double fwdFactor{ std::exp((mkt.r - mkt.q) * T) };

    // Convert plain moneyness (K/S0) into log-moneyness log(K/F)
    std::vector<double> logFM{};
    logFM.reserve(mny.size());
    std::transform(mny.begin(), mny.end(),
        std::back_inserter(logFM), [fwdFactor](double mny)
        {
            return std::log(mny / fwdFactor);
        });

    // Convert implied volatility into total variance
    std::vector<double> tVar{};
    tVar.reserve(vol.size());
    std::transform(vol.begin(), vol.end(),
        std::back_inserter(tVar), [T](double sigma)
        {
            return sigma * sigma * T;
        });

    return SliceData(T, mny, logFM, vol, tVar, mkt);
}

SliceData SliceData::fromModelData(const std::vector<double>& logFM,
    const std::vector<double>& wT,
    const MarketData& mkt,
    double T)
{
    // Compute forward factor: (r - q) * T
    double fwdFactor{ (mkt.r - mkt.q) * T };

    // Convert log-moneyness log(K/F) into plain moneyness K/S
    std::vector<double> mny{};
    mny.reserve(logFM.size());
    std::transform(logFM.begin(), logFM.end(),
        std::back_inserter(mny), [fwdFactor](double mny)
        {
            return std::exp(mny + fwdFactor);
        });

    // Convert total variance into implied volatility
    std::vector<double> vol{};
    vol.reserve(wT.size());
    std::transform(wT.begin(), wT.end(),
        std::back_inserter(vol), [T](double var)
        {
            return std::sqrt(var / T);
        });

    return SliceData(T, mny, logFM, vol, wT, mkt);
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
    double r{ marketData_.r };
    double q{ marketData_.q };
    double S{ marketData_.S };
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