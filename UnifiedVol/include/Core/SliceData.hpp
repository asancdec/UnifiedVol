/**
* SliceData.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 09/09/2025
*/

#ifndef SLICE_DATA_HPP
#define SLICE_DATA_HPP

#include "Core/MarketData.hpp"

#include <vector>

// Class to hold volatility surface data for one maturity
class SliceData
{   
private:

    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------

    double T_;                      // slice maturity
    std::vector<double> mny_;            // plain moneyness (K/S)
    std::vector<double> logFM_;          // log-forward moneyness log(K/F)
    std::vector<double> vol_;            // implied volatilities 
    std::vector<double> wT_;             // total variance (vol²T) 


    // Custom constructor using plain moneyness and implied volatility data
    explicit SliceData (double T,
                        const std::vector<double>& mny, 
                        const std::vector<double>& logFM,
                        const std::vector<double>& vol, 
                        const std::vector<double>& wT);
public:

    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------   
    SliceData() = delete;

    // Initialize class using plain moneyness and implied volatility data
    static SliceData fromMarketData(const std::vector<double>& mny,
        const std::vector<double>& vol,
        const MarketData& mkt,
        double T);

    // Initialize class using log-forward moneyness and total implied variance
    static SliceData fromModelData(const std::vector<double>& logFM,
        const std::vector<double>& wT,
        const MarketData& mkt,
        double T);

    //--------------------------------------------------------------------------
    // Math functions
    //--------------------------------------------------------------------------

    // Determine minimum total variance
    const double minWT() const noexcept;

    // Determine maximum total variance
    const double maxWT() const noexcept;

    // Determine minimum log-forward moneyness
    const double minLogFM() const noexcept;

    // Determine maximum log-forward moneyness
    const double maxLogFM() const noexcept;

    // Return ATM total variance
    const double atmWT() const noexcept;

    //--------------------------------------------------------------------------
    // Utilities
    //--------------------------------------------------------------------------
    
    // Print implied volatility slice on console
    void printImpVol() const noexcept;

    // Print total variance slice on console
    void printTotVar() const noexcept;

    // Getters
    const double T() const noexcept;
    const std::vector<double>& mny() const noexcept;
    const std::vector<double>& logFM() const noexcept;
    const std::vector<double>& vol() const noexcept;
    const std::vector<double>& wT() const noexcept;

};

#endif // SLICE_DATA_HPP