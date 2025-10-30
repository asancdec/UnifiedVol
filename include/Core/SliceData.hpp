/**
* SliceData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/MarketData.hpp"

#include <vector>

// Class to hold volatility surface data for one maturity
class SliceData
{   
private:

    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------

    double T_;                           // slice maturity
    double F_;                           // forward price
    std::vector<double> mny_;            // plain moneyness (K/S)
    std::vector<double> logFM_;          // log-forward moneyness log(K/F)
    std::vector<double> vol_;            // volatilities 
    std::vector<double> wT_;             // total variance (vol²T) 
    std::vector<double> K_;              // strikes vector
    std::vector<double> callBS_;         // value of European call options
    MarketData mktData_;                  // market data struct

public:

    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------   
    SliceData() = delete;
    explicit SliceData(double T,
        const std::vector<double>& mny,
        const std::vector<double>& vol,
        const MarketData& mktData);

    //--------------------------------------------------------------------------
    // Math functions
    //--------------------------------------------------------------------------

    // Black-Scholes pricing
    static double blackScholes(double T,
        double r,
        double q,
        double vol,
        double S,
        double K, 
        bool isCall) noexcept;

    // Normal cumulative density function
    static double normalCDF(double x) noexcept;

    // Determine minimum total variance
    double minWT() const noexcept;

    // Determine maximum total variance
    double maxWT() const noexcept;

    // Determine minimum log-forward moneyness
    double minLogFM() const noexcept;

    // Determine maximum log-forward moneyness
    double maxLogFM() const noexcept;

    // Return ATM total variance
    double atmWT() const noexcept;

    //--------------------------------------------------------------------------
    // Getters
    //--------------------------------------------------------------------------
 
    double T() const noexcept;
    double F() const noexcept;
    const std::vector<double>& mny() const noexcept;
    const std::vector<double>& logFM() const noexcept;
    const std::vector<double>& vol() const noexcept;
    const std::vector<double>& wT() const noexcept;
    const std::vector<double>& K() const noexcept;
    const std::vector<double>& callBS() const noexcept;


    //--------------------------------------------------------------------------
    // Setters
    //--------------------------------------------------------------------------

    void setWT(const std::vector<double>& wT); 
};

