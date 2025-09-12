/**
* SliceData.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/9/2025
*/

#include "Core/SliceData.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>


SliceData::SliceData(const std::vector<double>& m,
    const std::vector<double>& logFM,
    const std::vector<double>& vol, 
    const std::vector<double>& wT)
    : m(m), logFM(logFM), vol(vol), wT(wT)
{}

SliceData SliceData::FromMarketData(const std::vector<double>& m, 
    const std::vector<double>& vol,
    const MarketData& mkt, double T)
{   
    // Compute forward factor F/S0
    double fwdFactor{ std::exp((mkt.r - mkt.q) * T) };

    // Convert plain moneyness (K/S0) into log-moneyness log(K/F)
    std::vector<double> logFM{};
    logFM.reserve(m.size());
    std::transform(m.begin(), m.end(),
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

    return SliceData(m, logFM, vol, tVar);
}


SliceData SliceData::FromModelData(const std::vector<double>& logFM,
    const std::vector<double>& wT,
    const MarketData& mkt,
    double T)
{
    // Compute forward factor: (r - q) * T
    double fwdFactor{ (mkt.r - mkt.q) * T };

    // Convert log-moneyness log(K/F) into plain moneyness K/S
    std::vector<double> m{};
    m.reserve(logFM.size());
    std::transform(logFM.begin(), logFM.end(),
        std::back_inserter(m), [fwdFactor](double mny)
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

    return SliceData(m, logFM, vol, wT);
}


const double SliceData::minWT() const noexcept
{
    return *std::min_element(wT.begin(), wT.end());
}

const double SliceData::maxWT() const noexcept
{
    return *std::max_element(wT.begin(), wT.end());
}

const double SliceData::minLogFM() const noexcept
{
    return *std::min_element(logFM.begin(), logFM.end());
}

const double SliceData::maxLogFM() const noexcept
{
    return *std::max_element(logFM.begin(), logFM.end());
}

void SliceData::printImpVol() const
{
    for (const auto& v : vol)
    {
        std::cout << std::fixed << std::setprecision(4) << v << "\t";
    }
    std::cout << "\n";
}

void SliceData::printTotVar() const
{
    for (const auto& v : wT)
    {
        std::cout << std::fixed << std::setprecision(4) << v << "\t";
    }
    std::cout << "\n";
}