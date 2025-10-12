/**
* VolSurface.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Core/VolSurface.hpp"
#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <string>

using uv::ErrorCode;

VolSurface::VolSurface(std::vector<SliceData>&& slices,
    std::vector<double>&& maturities)
: slices_(std::move(slices)), maturities_(std::move(maturities)) {}

VolSurface VolSurface::fromMarketData(const std::vector<double>& mny,
    const std::vector<std::vector<double>>& vols,
    const std::vector<double>& maturities,
    const MarketData& mkt)
{
    // --- Size checks ---
    UV_REQUIRE(vols.size() == maturities.size(),
        ErrorCode::InvalidArgument,
        "vols.size() != maturities.size()");
    for (size_t i = 0; i < vols.size(); ++i) 
    {
        UV_REQUIRE(vols[i].size() == mny.size(),
            ErrorCode::InvalidArgument,
            "vols[" + std::to_string(i) + "].size() != mny.size()");
    }

    // --- Build slices ---
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
    // --- Size checks ---
    UV_REQUIRE(logFM.size() == maturities.size(),
        ErrorCode::InvalidArgument,
        "logFM.size() != maturities.size()");
    UV_REQUIRE(wT.size() == maturities.size(),
        ErrorCode::InvalidArgument,
        "wT.size() != maturities.size()");
    for (size_t i = 0; i < maturities.size(); ++i) {
        UV_REQUIRE(logFM[i].size() == wT[i].size(),
            ErrorCode::InvalidArgument,
            "logFM[" + std::to_string(i) + "].size() != wT[" + std::to_string(i) + "].size()");
    }

    // --- Build slices ---
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
    std::ostringstream oss;
    oss << '\n';
    oss << "T\\%S\t";

    // Header row (moneyness)
    for (const auto& m : slices_[0].mny())
        oss << std::fixed << std::setprecision(2) << m << '\t';
    oss << '\n';

    // Each maturity row
    for (size_t i = 0; i < slices_.size(); ++i)
    {
        oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

        for (const auto& v : slices_[i].vol())
            oss << std::fixed << std::setprecision(4) << v << '\t';

        oss << '\n';
    }

    UV_INFO(oss.str()); 
}

void VolSurface::printTotVar() const noexcept
{
    std::ostringstream oss;
    oss << '\n';
    oss << "T\\k\t";

    // Header row (moneyness)
    for (const auto& m : slices_[0].mny())
        oss << std::fixed << std::setprecision(2) << m << '\t';
    oss << '\n';

    // Each maturity row
    for (size_t i = 0; i < slices_.size(); ++i)
    {
        oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

        for (const auto& v : slices_[i].wT())
            oss << std::fixed << std::setprecision(4) << v << '\t';

        oss << '\n';
    }

    UV_INFO(oss.str()); 
}

void VolSurface::printBSCall() const noexcept
{
    std::ostringstream oss;
    oss << '\n';
    oss << "T\\k\t";

    // Header row (moneyness)
    for (const auto& m : slices_[0].mny())
        oss << std::fixed << std::setprecision(2) << m << '\t';
    oss << '\n';

    // Each maturity row
    for (size_t i = 0; i < slices_.size(); ++i)
    {
        oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

        for (const auto& v : slices_[i].callBS())
            oss << std::fixed << std::setprecision(2) << v << '\t';

        oss << '\n';
    }
    UV_INFO(oss.str());
}

std::size_t VolSurface::numStrikes() const
{   
    // Extract the number of strikes in the first slice
    const std::size_t n{ slices_.front().logFM().size() };

    // Check that every slice has the same number of strikes
    for (std::size_t i = 1; i < slices_.size(); ++i)
    {
        const std::size_t ni{ slices_[i].logFM().size() };
        UV_REQUIRE(ni == n, ErrorCode::InvalidArgument,
            "numStrikes(): inconsistent k-grid length — slice 0 has " +
            std::to_string(n) + " strikes, slice " + std::to_string(i) +
            " has " + std::to_string(ni));
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
