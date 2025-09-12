/**
* VolSurface.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/


#ifndef VOLSURFACE_HPP
#define VOLSURFACE_HPP

#include <vector>
#include <stdexcept>

class VolSurface
{
public:

    // Public attributes
    std::vector<double> maturities; // Annualized maturities
    std::vector<double> moneyness;    // Moneyness (K/S)
    std::vector<std::vector<double>> vols; // 2D volatility grid: vols[maturity][moneyness]
    double r;  // Continously compounded annualized risk-free rate
    double q;  // Continously compounded annualized dividend yield

    // Delete default constructor
    VolSurface() = delete;

    // Custom constructor
    VolSurface(const std::vector<double>& moneynessInput,
               const std::vector<double>& maturitiesInput,
               const std::vector<std::vector<double>>& volsInput,
               double rInput = 0,
               double qInput = 0);

    // Function to get volatility surface dimensions
    std::pair<size_t, size_t> getDimensions() const;

    // Calculate log forward moneyness -> log(K/S)
    std::vector<std::vector<double>> k() const;

    // Calculate implied variance time
    std::vector<std::vector<double>> wMkt() const;

    // Find ATM index and the indexes before and after
    std::tuple<size_t, size_t, size_t> findATMIndices() const;

    // Print volatility surface on the console
    void printConsole() const;
};


#endif // VOLSURFACE_HPP