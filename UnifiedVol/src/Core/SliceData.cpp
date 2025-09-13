/**
* SliceData.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/9/2025
*/

#include "Core/SliceData.hpp"

#include <cmath>
#include <numbers>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>


SliceData::SliceData(double T,
    const std::vector<double>& mny,
    const std::vector<double>& logFM,
    const std::vector<double>& vol, 
    const std::vector<double>& wT)
    : T_(T), mny_(mny), logFM_(logFM), vol_(vol), wT_(wT)
{}

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

    return SliceData(T, mny, logFM, vol, tVar);
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

    return SliceData(T, mny, logFM, vol, wT);
}

const double SliceData::minWT() const noexcept
{
    return *std::min_element(wT_.begin(), wT_.end());
}

const double SliceData::maxWT() const noexcept
{
    return *std::max_element(wT_.begin(), wT_.end());
}

const double SliceData::minLogFM() const noexcept
{
    return *std::min_element(logFM_.begin(), logFM_.end());
}

const double SliceData::maxLogFM() const noexcept
{
    return *std::max_element(logFM_.begin(), logFM_.end());
}

const double SliceData::atmWT() const noexcept
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

std::vector<double> SliceData::vega() const noexcept
{
    const double sqrtT{ std::sqrt(T_) };
    std::vector<double> v(logFM_.size());

    for (size_t i = 0; i < logFM_.size(); ++i)
    {
        const double d1{ (-logFM_[i] + 0.5 * vol_[i] * vol_[i] * T_) / (vol_[i] * sqrtT) };
        const double nd1{ std::exp(-0.5 * d1 * d1) / std::sqrt(2.0 * std::numbers::pi) };
        v[i] = nd1 * sqrtT;
    }
    return v;
}


void SliceData::printImpVol() const noexcept
{
    for (const auto& v : vol_)
    {
        std::cout << std::fixed << std::setprecision(4) << v << "\t";
    }
    std::cout << "\n";
}

void SliceData::printTotVar() const noexcept
{   
    for (const auto& v : wT_)
    {
        std::cout << std::fixed << std::setprecision(4) << v << "\t";
    }
    std::cout << "\n";
}

const double SliceData::T() const noexcept
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


