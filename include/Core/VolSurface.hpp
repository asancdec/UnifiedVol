/**
* VolSurface.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/SliceData.hpp"
#include "Core/MarketData.hpp"

#include <vector>

namespace uv
{
    class VolSurface
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        ::std::vector<SliceData> slices_;         // 2D volatility grid: vols[maturity][maturity slice]
        ::std::vector<double>    maturities_;     // Vector with different volatility surface maturities

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        VolSurface() = delete;
        explicit VolSurface(const ::std::vector<double>& mny,
            const ::std::vector<::std::vector<double>>& vols,
            const ::std::vector<double>& maturities,
            const MarketData& mktData);

        //--------------------------------------------------------------------------
        // Utilities
        //--------------------------------------------------------------------------

        // Print volatility surface on the console
        void printVol() const noexcept;

        // Print total variance surface on the console
        void printTotVar() const noexcept;

        // Print Black-Scholes Call price surface
        void printBSCall() const noexcept;

        // Return number of strikes
        ::std::size_t numStrikes() const;

        // Getters
        ::std::vector<SliceData>& slices() noexcept;
        const ::std::vector<double>& maturities() const noexcept;
        size_t numSlices() const noexcept;
    };
}
