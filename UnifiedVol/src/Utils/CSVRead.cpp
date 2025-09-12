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


VolSurface readVolSurface(const std::string& filename)
{
    // Open the CSV file for reading
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << "\n";
        // Could consider throwing an exception here instead of continuing
    }

    // Container to store the raw CSV data
    std::vector<std::vector<double>> data{};
    std::string line{};

    // Read file line by line
    while (std::getline(file, line))
    {
        std::stringstream ss(line); // Split line by commas
        std::string cell{};
        std::vector<double> row{}; // Row of doubles

        // Read each cell in the line
        while (std::getline(ss, cell, ','))
        {
            try
            {
                row.push_back(std::stod(cell)); // Convert string to double
            }
            catch (...)
            {
                // Ignore any non-numeric cells (e.g., empty or header)
            }
        }

        // Only add non-empty rows to the data
        if (!row.empty())
        {
            data.push_back(row);
        }
    }

    // Extract moneyness (K/S) from the first row
    // Skip the first element (row[0]) because it may be an empty header cell
    std::vector<double> moneyness(data[0].begin(), data[0].end());

    // Prepare vector for maturities
    std::vector<double> maturities{};
    maturities.reserve(data.size() - 1); 

    // Prepare matrix for implied volatilities
    std::vector<std::vector<double>> vols{};
    vols.reserve(data.size() - 1);

    // Loop through each row (skip first row because it's the header)
    for (size_t i = 1; i < data.size(); ++i)
    {
        // First column in the row is the maturity
        maturities.push_back(data[i][0]);

        // Remaining columns are the implied volatilities
        std::vector<double> rowVols(data[i].begin() + 1, data[i].end());
        vols.push_back(rowVols);
    }

    // Construct and return a VolSurface object with the extracted data
    return VolSurface(moneyness, maturities, vols);
}