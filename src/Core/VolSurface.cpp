/**
* VolSurface.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Core/VolSurface.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <string>

namespace uv::core
{
    VolSurface::VolSurface(const std::vector<double>& mny,
        const std::vector<std::vector<double>>& vols,
        const std::vector<double>& tenors,
        const MarketData& mktData)
        : tenors_(tenors),
        numTenors_(tenors.size())
    {   
        // ---------- Sanity checks ------
        
        // Number of volatility slices (one per maturity)
        const std::size_t dim0{ vols.size()};

        // Validate that the volatility matrix provides one slice per maturity
        UV_REQUIRE(
            dim0 == numTenors_,
            ErrorCode::InvalidArgument,
            "VolSurface: number of vol slices (" + std::to_string(dim0) +
            ") does not match number of tenors (" + std::to_string(numTenors_) + ")"
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
        slices_.reserve(numTenors_);
        for (size_t i = 0; i < numTenors_; ++i)
        {
            slices_.emplace_back
            (
                SliceData(tenors_[i],
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

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

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

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

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

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : slices_[i].callBS())
                oss << std::fixed << std::setprecision(3) << v << '\t';

            oss << '\n';
        }
        UV_INFO(oss.str());
    }

    std::vector<std::vector<double>> VolSurface::totVarMatrix() const noexcept
    {
        // Allocate matrix:
        // outer dimension = tenors (rows)
        // inner dimension = strikes (columns)        
        std::vector<std::vector<double>> matrix(
            numTenors_, std::vector<double>(numStrikes_)
        );

        // Populate total variance slices
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            // Extract total variance slice at tenor i
            matrix[i] = slices_[i].wT();
        }

        return matrix;
    }

    std::vector<SliceData>& VolSurface::slices() noexcept
    {
        return slices_;
    }

    const std::vector<double>& VolSurface::tenors() const noexcept
    {
        return tenors_;
    }

    std::size_t VolSurface::numTenors() const noexcept
    {
        return numTenors_;
    }

    std::size_t VolSurface::numStrikes() const noexcept
    {
        return numStrikes_;
    }

}
