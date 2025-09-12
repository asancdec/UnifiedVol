/**
* VolSurface.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/21/2025
*/

#include "Core/VolSurface.hpp"
#include <iomanip>
#include <iostream>
#include <algorithm>

VolSurface::VolSurface(std::vector<SliceData>&& slices,
    std::vector<double>&& mats)
: slices(std::move(slices)), mats(std::move(mats)) {}


VolSurface VolSurface::FromMarketData(const std::vector<double>& m,
    const std::vector<std::vector<double>>& vols,
    const std::vector<double>& mats,
    const MarketData& mkt)
{
    std::vector<SliceData> slices;
    slices.reserve(mats.size());

    for (size_t i = 0; i < mats.size(); ++i) 
    {
        slices.push_back
        (
            SliceData::FromMarketData(m, vols[i], mkt, mats[i])
        );
    }
    return VolSurface(std::move(slices), std::vector<double>(mats));
}


VolSurface VolSurface::FromModelData(const std::vector<std::vector<double>>& logFM,
    const std::vector<std::vector<double>>& wT,
    const std::vector<double>& mats,
    const MarketData& mkt)
{
    std::vector<SliceData> slices;
    slices.reserve(mats.size());

    for (size_t i = 0; i < mats.size(); ++i)
    {
        slices.push_back
        (
            SliceData::FromModelData(logFM[i], wT[i], mkt, mats[i])
        );
    }
    return VolSurface(std::move(slices), std::vector<double>(mats));
}


void VolSurface::printImpVol() const
{
    // Print header row
    std::cout << "T\\%S\t";

    // Print moneyness K/S from the first slice, 2 decimals
    for (const auto& m : slices[0].m)
    {
        std::cout << std::fixed << std::setprecision(2) << m << "\t";
    }
    std::cout << "\n";

    // Print implied volatility slices
    for (size_t i = 0; i < slices.size(); ++i)
    {
        std::cout << std::fixed << std::setprecision(2) << mats[i] << "\t";
        slices[i].printImpVol();
    }
}

void VolSurface::printTotVar() const
{
    // Print header row
    std::cout << "T\\k\t";

    // Print moneyness K/S from the first slice, 2 decimals
    for (const auto& m : slices[0].m)
    {
        std::cout << std::fixed << std::setprecision(2) << m << "\t";
    }
    std::cout << "\n";

    // Print total variance volatility slices
    for (size_t i = 0; i < slices.size(); ++i)
    {
        std::cout << std::fixed << std::setprecision(2) << mats[i] << "\t";
        slices[i].printTotVar();
    }

    // Print total variance volatility slices
    for (size_t i = 0; i < slices.size(); ++i)
    {
        std::cout << "\n";
        std::cout<< slices[i].maxLogFM();
    }
}

//
//// Find ATM index and the indexes before and after
//std::tuple<size_t, size_t, size_t> VolSurface::findATMIndices() const
//{
//    auto it = std::min_element(
//        this->moneyness.begin(), this->moneyness.end(),
//        [](double k1, double k2) { return std::abs(k1 - 1.0) < std::abs(k2 - 1.0); }
//    );
//
//    size_t atmIndex = std::distance(this->moneyness.begin(), it);
//    size_t prevIndex = (atmIndex > 0) ? atmIndex - 1 : atmIndex;
//    size_t nextIndex = (atmIndex + 1 < this->moneyness.size()) ? atmIndex + 1 : atmIndex;
//
//    return { atmIndex, prevIndex, nextIndex };
//}