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
    VolSurface::VolSurface(const std::vector<double>& mny,
        const std::vector<std::vector<double>>& vols,
        const std::vector<double>& maturities,
        const MarketData& mktData)
        : maturities_(maturities),
        numMaturities_(maturities_.size())
    {   
        // ---------- Sanity checks ------
        
        // Number of volatility slices (one per maturity)
        const std::size_t dim0{ vols.size()};

        // Validate that the volatility matrix provides one slice per maturity
        UV_REQUIRE(
            dim0 == numMaturities_,
            ErrorCode::InvalidArgument,
            "VolSurface: number of vol slices (" + std::to_string(dim0) +
            ") does not match number of maturities (" + std::to_string(numMaturities_) + ")"
        );

        // Number of strikes in the strike grid (must be consistent across slices)
        const std::size_t dim1{ mny.size() };

        // Validate that each volatility slice has the same number of strikes
        for (std::size_t i = 0; i < dim0; ++i)
        {
            const std::size_t t{ vols[i].size() };
            UV_REQUIRE(
                t == dim1,
                ErrorCode::InvalidArgument,
                "VolSurface: inconsistent strike dimension — expected " +
                std::to_string(t) + " strikes, but slice " + std::to_string(i) +
                " has " + std::to_string(dim1)
            );
        }

        // ---------- Initialize member variables ------
        numStrikes_ = dim1;
        slices_.reserve(numMaturities_);
        for (size_t i = 0; i < numMaturities_; ++i)
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
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\%S\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each maturity row
        for (size_t i = 0; i < numMaturities_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].vol())
                oss << std::fixed << std::setprecision(5) << v << '\t';

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
        for (size_t i = 0; i < numMaturities_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].wT())
                oss << std::fixed << std::setprecision(5) << v << '\t';

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
        for (size_t i = 0; i < numMaturities_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << maturities_[i] << '\t';

            for (const auto& v : slices_[i].callBS())
                oss << std::fixed << std::setprecision(3) << v << '\t';

            oss << '\n';
        }
        UV_INFO(oss.str());
    }

    std::vector<std::vector<double>> VolSurface::totVarMatrix() const noexcept
    {
        // Allocate matrix:
        // outer dimension = maturities (rows)
        // inner dimension = strikes (columns)        
        std::vector<std::vector<double>> matrix(
            numMaturities_, std::vector<double>(numStrikes_)
        );

        // Populate total variance slices
        for (std::size_t i = 0; i < numMaturities_; ++i)
        {
            // Extract total variance slice at maturity i
            matrix[i] = slices_[i].wT();
        }

        return matrix;
    }

    std::vector<SliceData>& VolSurface::slices() noexcept
    {
        return slices_;
    }

    const std::vector<double>& VolSurface::maturities() const noexcept
    {
        return maturities_;
    }

    std::size_t VolSurface::numMaturities() const noexcept
    {
        return numMaturities_;
    }

    std::size_t VolSurface::numStrikes() const noexcept
    {
        return numStrikes_;
    }

}
