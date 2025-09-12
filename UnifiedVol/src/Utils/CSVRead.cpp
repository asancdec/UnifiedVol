/**
* CSVRead.cpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/

#include "Utils/CSVRead.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>


VolSurface readVolSurface(const std::string& filename, const MarketData& mkt)
{
    // Open the CSV file for reading
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // Container to store the raw CSV data
    std::vector<std::vector<double>> data{};
    std::string line{};

    // Read file line by line
    while (std::getline(file, line)) 
    {
        std::stringstream ss(line);
        std::string cell{};
        std::vector<double> row{};
        while (std::getline(ss, cell, ','))
        {
            try 
            {
                row.push_back(std::stod(cell));
            }
            catch (...) 
            {
                // Ignore non-numeric cells (e.g., empty or header)
            }
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }

    if (data.empty()) {
        throw std::runtime_error("CSV file is empty or invalid: " + filename);
    }

    // Extract moneyness (K/S) from the first row
    std::vector<double> mny(data[0].begin(), data[0].end());

    // Prepare maturities and vol matrix
    std::vector<double> maturities{};
    std::vector<std::vector<double>> vols{};
    maturities.reserve(data.size() - 1);
    vols.reserve(data.size() - 1);

    for (size_t i = 1; i < data.size(); ++i) 
    {
        maturities.push_back(data[i][0]); // First column is maturity
        vols.push_back(std::vector<double>(data[i].begin() + 1, data[i].end())); // Remaining are vols
    }

    // Construct VolSurface using the new FromMarketData method
    return VolSurface::fromMarketData(mny, vols, maturities, mkt);
}