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

    // Custom constructor using plain moneyness and implied volatility data
    explicit SliceData (const std::vector<double>& m, 
                        const std::vector<double>& logFM,
                        const std::vector<double>& vol, 
                        const std::vector<double>& wT);
     
public:
    std::vector<double> m;            // plain moneyness (K/S)
    std::vector<double> logFM;        // log-forward moneyness log(K/F)
    std::vector<double> vol;          // implied volatilities 
    std::vector<double> wT;           // total variance (vol²T) 

    // Initialization
   
    SliceData() = delete;

    // Initialize class using plain moneyness and implied volatility data
    static SliceData FromMarketData(const std::vector<double>& m,
        const std::vector<double>& vol,
        const MarketData& mkt,
        double T);

    // Initialize class using log-forward moneyness and total implied variance
    static SliceData FromModelData(const std::vector<double>& logFM,
        const std::vector<double>& wT,
        const MarketData& mkt,
        double T);

    // Functionalities

    // Determine minimum total variance
    const double minWT() const noexcept;

    // Determine maximum total variance
    const double maxWT() const noexcept;

    // Determine minimum log-forward moneyness
    const double minLogFM() const noexcept;

    // Determine maximum log-forward moneyness
    const double maxLogFM() const noexcept;

    // Utilities
    
    // Print implied volatility slice on console
    void printImpVol() const;
    // Print total variance slice on console
    void printTotVar() const;
};

#endif // SLICE_DATA_HPP