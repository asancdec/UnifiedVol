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

namespace uv 
{
    VolSurface::VolSurface(const ::std::vector<double>& mny,
        const ::std::vector<::std::vector<double>>& vols,
        const ::std::vector<double>& maturities,
        const MarketData& mktData)
        : maturities_(maturities)
    {
        // Build slices
        slices_.reserve(maturities.size());
        for (size_t i = 0; i < maturities.size(); ++i)
        {
            slices_.emplace_back
            (
                SliceData(maturities_[i],
                    mny,
                    vols[i],
                    mktData)
            );
        }
    }

    void VolSurface::printVol() const noexcept
    {
        ::std::ostringstream oss;
        oss << '\n';
        oss << "T\\%S\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << ::std::fixed << ::std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each maturity row
        for (size_t i = 0; i < slices_.size(); ++i)
        {
            oss << ::std::fixed << ::std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].vol())
                oss << ::std::fixed << ::std::setprecision(5) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printTotVar() const noexcept
    {
        ::std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << ::std::fixed << ::std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each maturity row
        for (size_t i = 0; i < slices_.size(); ++i)
        {
            oss << ::std::fixed << ::std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].wT())
                oss << ::std::fixed << ::std::setprecision(4) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printBSCall() const noexcept
    {
        ::std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << ::std::fixed << ::std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each maturity row
        for (size_t i = 0; i < slices_.size(); ++i)
        {
            oss << ::std::fixed << ::std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].callBS())
                oss << ::std::fixed << ::std::setprecision(2) << v << '\t';

            oss << '\n';
        }
        UV_INFO(oss.str());
    }

    ::std::size_t VolSurface::numStrikes() const
    {
        // Extract the number of strikes in the first slice
        const ::std::size_t n{ slices_.front().logFM().size() };

        // Check that every slice has the same number of strikes
        for (::std::size_t i = 1; i < slices_.size(); ++i)
        {
            const ::std::size_t ni{ slices_[i].logFM().size() };
            UV_REQUIRE(ni == n, ErrorCode::InvalidArgument,
                "numStrikes(): inconsistent k-grid length — slice 0 has " +
                ::std::to_string(n) + " strikes, slice " + ::std::to_string(i) +
                " has " + ::std::to_string(ni));
        }
        return n;
    }

    ::std::vector<SliceData>& VolSurface::slices() noexcept
    {
        return slices_;
    }

    const ::std::vector<double>& VolSurface::maturities() const noexcept
    {
        return maturities_;
    }

    size_t VolSurface::numSlices() const noexcept
    {
        return slices_.size();
    }
}
