/**
* CSVRead.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/VolSurface.hpp"

#include <vector>
#include <string>

namespace uv::utils
{
/**
 * @brief Reads a CSV file containing a volatility surface and returns a VolSurface object.
 *
 * The CSV file is expected to have the following structure:
 * - The first row contains moneyness (K/S).
 * - The first column of each subsequent row contains tenors.
 * - The remaining cells contain implied volatilities for each strike-maturity pair.
 *
 * @param filename Path to the CSV file.
 * @return VolSurface Object containing the strikes, tenors, and implied volatilities.
 */
VolSurface readVolSurface(const std::string& filename, const MarketData& mktData);
}

