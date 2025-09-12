/**
* VolSurface.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/


#include "DataStructs/VolSurface.hpp"
#include <iomanip>
#include <iostream>
#include <algorithm>

// Custom constructor
VolSurface::VolSurface(const std::vector<double>& moneynessInput,
    const std::vector<double>& maturitiesInput,
    const std::vector<std::vector<double>>& volsInput,
    double rInput,
    double qInput)
    : moneyness(moneynessInput), maturities(maturitiesInput), vols(volsInput), r(rInput), q(qInput)
{
    if (vols.size() != maturities.size())
        throw std::invalid_argument("Vol rows must match maturities size.");
    for (const auto& row : vols)
        if (row.size() != moneyness.size())
            throw std::invalid_argument("Vol columns must match moneyness size.");
}

std::pair<size_t, size_t> VolSurface::getDimensions() const
{
    return { this->maturities.size(), this->moneyness.size() };
}


// Method to calculate log-forward moneyness from the Implied Volatility surface
std::vector<std::vector<double>> VolSurface::k() const
{
    // Initialize vector to store log-forward moneyness calculations
    std::vector<std::vector<double>> k{};
    k.reserve(this->maturities.size());

    // Iterate through each maturity (row)
    for (const auto& T : this->maturities)
    {
        // Precompute forward 
        double fwd{ std::exp((this->r - this->q) * T) };

        std::vector<double> row{};
        row.reserve(this->moneyness.size());

        // Transform moneyness (K/S) into log-forward moneyness (log(K/F))
        std::transform(this->moneyness.begin(), this->moneyness.end(), std::back_inserter(row),
            [fwd](double moneyness) {return std::log(moneyness / fwd);  });

        k.push_back(std::move(row));
    }

    return k;
}


// Method to calculate implied time variance from the Implied Volatility surface
std::vector<std::vector<double>> VolSurface::wMkt() const
{

    // Initialize vector to store results
    std::vector<std::vector<double>> wMkt{};
    wMkt.reserve(this->maturities.size());

    // Iterate through each maturity (row)
    for (size_t i = 0; i < this->maturities.size(); ++i)
    {
        double T{ this->maturities[i] };

        std::vector<double> row{};
        row.reserve(this->vols[i].size());

        // Transform implied volatility into variance time
        std::transform(this->vols[i].begin(), this->vols[i].end(), std::back_inserter(row),
            [T](double sigma) { return sigma * sigma * T; });

        wMkt.push_back(std::move(row));
    }

    return wMkt;
}


void VolSurface::printConsole() const
{
    // Print volatility surface dimensions
    const auto [numRows, numCols] = getDimensions();
    std::cout << "VolSurface dimensions: " << numRows << " x " << numCols << "\n";

    // Print title row
    std::cout << "T\\K/S\t";

    // Print moneyness (K/S)
    for (const auto& s : this->moneyness)
        std::cout << std::fixed << std::setprecision(2) << s << "\t";
    std::cout << "\n";

    // Iterate from the second to the last row
    for (size_t i = 0; i < numRows; ++i)
    {
        // Print maturity with 2 decimals
        std::cout << std::fixed << std::setprecision(2) << this->maturities[i] << "\t";

        // Print vols with 4 decimals
        for (const auto& v : this->vols[i])
            std::cout << std::fixed << std::setprecision(4) << v << "\t";
        std::cout << "\n";
    }
}


// Find ATM index and the indexes before and after
std::tuple<size_t, size_t, size_t> VolSurface::findATMIndices() const
{
    auto it = std::min_element(
        this->moneyness.begin(), this->moneyness.end(),
        [](double k1, double k2) { return std::abs(k1 - 1.0) < std::abs(k2 - 1.0); }
    );

    size_t atmIndex = std::distance(this->moneyness.begin(), it);
    size_t prevIndex = (atmIndex > 0) ? atmIndex - 1 : atmIndex;
    size_t nextIndex = (atmIndex + 1 < this->moneyness.size()) ? atmIndex + 1 : atmIndex;

    return { atmIndex, prevIndex, nextIndex };
}