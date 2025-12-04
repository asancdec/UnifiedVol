/**
* VolSurface.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/SliceData.hpp"
#include "Core/MarketData.hpp"

#include <vector>
#include <cstddef>

namespace uv
{
    class VolSurface
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        std::vector<SliceData> slices_;         // 2D volatility grid: vols[maturity][maturity slice]
        std::vector<double>    maturities_;     // Vector with different volatility surface maturities
        std::size_t numMaturities_;             // Number of maturities
        std::size_t numStrikes_;                // Number of strikes

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        VolSurface() = delete;
        explicit VolSurface(const std::vector<double>& mny,
            const std::vector<std::vector<double>>& vols,
            const std::vector<double>& maturities,
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

        // Return total variance matrix
        std::vector<std::vector<double>> totVarMatrix() const noexcept;

        // Getters
        std::vector<SliceData>& slices() noexcept;
        const std::vector<double>& maturities() const noexcept;
        std::size_t numMaturities() const noexcept;
        std::size_t numStrikes() const noexcept;
    };
}
