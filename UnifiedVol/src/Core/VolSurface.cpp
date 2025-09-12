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
    std::vector<double>&& maturities)
: slices_(std::move(slices)), maturities_(std::move(maturities)) {}

VolSurface VolSurface::fromMarketData(const std::vector<double>& mny,
    const std::vector<std::vector<double>>& vols,
    const std::vector<double>& maturities,
    const MarketData& mkt)
{
    std::vector<SliceData> slices;
    slices.reserve(maturities.size());

    for (size_t i = 0; i < maturities.size(); ++i)
    {
        slices.push_back
        (
            SliceData::fromMarketData(mny, vols[i], mkt, maturities[i])
        );
    }
    return VolSurface(std::move(slices), std::vector<double>(maturities));
}

VolSurface VolSurface::fromModelData(const std::vector<std::vector<double>>& logFM,
    const std::vector<std::vector<double>>& wT,
    const std::vector<double>& maturities,
    const MarketData& mkt)
{
    std::vector<SliceData> slices;
    slices.reserve(maturities.size());

    for (size_t i = 0; i < maturities.size(); ++i)
    {
        slices.push_back
        (
            SliceData::fromModelData(logFM[i], wT[i], mkt, maturities[i])
        );
    }
    return VolSurface(std::move(slices), std::vector<double>(maturities));
}

void VolSurface::printImpVol() const
{
    // Print header row
    std::cout << "T\\%S\t";

    // Print moneyness K/S from the first slice, 2 decimals
    for (const auto& m : slices_[0].mny())
    {
        std::cout << std::fixed << std::setprecision(2) << m << "\t";
    }
    std::cout << "\n";

    // Print implied volatility slices
    for (size_t i = 0; i < slices_.size(); ++i)
    {
        std::cout << std::fixed << std::setprecision(2) << maturities_[i] << "\t";
        slices_[i].printImpVol();
    }
}

void VolSurface::printTotVar() const
{
    // Print header row
    std::cout << "T\\k\t";

    // Print moneyness K/S from the first slice, 2 decimals
    for (const auto& m : slices_[0].mny())
    {
        std::cout << std::fixed << std::setprecision(2) << m << "\t";
    }
    std::cout << "\n";

    // Print total variance slices
    for (size_t i = 0; i < slices_.size(); ++i)
    {
        std::cout << std::fixed << std::setprecision(2) << maturities_[i] << "\t";
        slices_[i].printTotVar();
    }
}

const std::vector<SliceData>& VolSurface::slices() const
{
    return slices_;
}

const std::vector<double>& VolSurface::maturities() const
{
    return maturities_;
}