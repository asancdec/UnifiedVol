/**
* VolSurface.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/

#ifndef VOLSURFACE_HPP
#define VOLSURFACE_HPP

#include "Core/SliceData.hpp"
#include "Core/MarketData.hpp"

#include <vector>

class VolSurface
{
private:

    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------

    std::vector<SliceData> slices_;         // 2D volatility grid: vols[maturity][maturity slice]
    std::vector<double>    maturities_;     // Vector with different volatility surface maturities

    // Custom constructor
    VolSurface(std::vector<SliceData>&& slices, std::vector<double>&& mats);

public:

    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    VolSurface() = delete;

    // Initialize class using plain moneyness K/S and implied volatility data
    static VolSurface fromMarketData(const std::vector<double>& mny,
        const std::vector<std::vector<double>>& vols,
        const std::vector<double>& mats,
        const MarketData& mkt);

    // Initialize class using log-forward moneyness log(K/F) and total variance
    static VolSurface fromModelData(const std::vector<std::vector<double>>& logFM,
        const std::vector<std::vector<double>>& wT,
        const std::vector<double>& maturities,
        const MarketData& mkt);

    //--------------------------------------------------------------------------
    // Utilities
    //--------------------------------------------------------------------------

    // Print volatility surface on the console
    void printVol() const noexcept;

    // Print total variance surface on the console
    void printTotVar() const noexcept;

    // Return number of strikes
    std::size_t numStrikes() const;

    // Getters
    std::vector<SliceData>& slices() noexcept;
    const std::vector<double>& maturities() const noexcept;
    size_t numSlices() const noexcept;
};


#endif // VOLSURFACE_HPP