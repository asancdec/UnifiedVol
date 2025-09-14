/**
* VolSurface.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/21/2025
*/

#include "Core/VolSurface.hpp"

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdexcept>


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

void VolSurface::printVol() const noexcept
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
        slices_[i].printVol();
    }
}

void VolSurface::printTotVar() const noexcept
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

std::size_t VolSurface::numStrikes() const
{   
    // Extract the number of strikes in the first slice
    const std::size_t n{ slices_.front().logFM().size() };

    // Check that every slice has the same number of strikes
    for (std::size_t i = 1; i < slices_.size(); ++i)
    {
        const std::size_t ni{ slices_[i].logFM().size() };
        if (ni != n) 
        {
            throw std::runtime_error
            (
                "numStrikes(): inconsistent k-grid length — slice 0 has " +
                std::to_string(n) + " strikes, slice " + std::to_string(i) +
                " has " + std::to_string(ni));
        }
    }
    return n;
}

std::vector<SliceData>& VolSurface::slices() noexcept
{
    return slices_;
}

const std::vector<double>& VolSurface::maturities() const noexcept
{
    return maturities_;
}

size_t VolSurface::numSlices() const noexcept
{
    return slices_.size();
}
