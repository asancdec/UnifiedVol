///**
// * CSVWrite.hpp
// * Author: Alvaro Sanchez de Carlos
// * Date: 08/18/2025
// */
//
//#ifndef CSV_WRITE_HPP
//#define CSV_WRITE_HPP
//
//#include "Core/VolSurface.hpp"
//#include <string>
//
// /**
//  * @brief Writes a VolSurface object to a CSV file.
//  *
//  * The CSV file will have the following structure:
//  * - The first row contains moneyness (K/S) headers.
//  * - The first column of each subsequent row contains maturities.
//  * - The remaining cells contain implied volatilities for each strike-maturity pair.
//  *
//  * @param volSurf  The VolSurface object to write.
//  * @param filename Path to the output CSV file.
//  */
//void writeVolSurface(const VolSurface& volSurf, 
//	const std::string& filename);
//
//#endif // CSVWRITE_HPP