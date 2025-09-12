/**
* CSVRead.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/17/2025
*/


#ifndef CSVREAD_HPP
#define CSVREAD_HPP

#include "DataStructs/VolSurface.hpp"
#include <vector>
#include <string>

/**
 * @brief Reads a CSV file containing a volatility surface and returns a VolSurface object.
 *
 * The CSV file is expected to have the following structure:
 * - The first row contains moneyness (K/S).
 * - The first column of each subsequent row contains maturities.
 * - The remaining cells contain implied volatilities for each strike-maturity pair.
 *
 * @param filename Path to the CSV file.
 * @return VolSurface<double> Object containing the strikes, maturities, and implied volatilities.
 */
VolSurface readVolSurface(const std::string& filename);



#endif // CSVREAD_HPP