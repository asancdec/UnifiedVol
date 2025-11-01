///**
// * CSVWrite.cpp
// * Author: Alvaro Sanchez de Carlos
// * Date: 08/18/2025
// */
//
//#include "Utils/CSVWrite.hpp"
//#include <fstream>
//#include <iostream>
//#include <iomanip>
//
//void writeVolSurface(const VolSurface& volSurf,
//    const ::std::string& filename)
//{
//    // Open the CSV file for writing
//    ::std::ofstream file(filename);
//    if (!file.is_open())
//    {
//        ::std::cerr << "Error: Could not open file " << filename << " for writing\n";
//        return; // Could throw exception instead
//    }
//
//    // Write header row (moneyness)
//    file << "Maturity/Moneyness";
//    for (const auto& k : volSurf.moneyness)
//    {
//        file << "," << ::std::setprecision(6) << k;
//    }
//    file << "\n";
//
//    // Write each row: maturity + implied volatilities
//    for (size_t i = 0; i < volSurf.maturities.size(); ++i)
//    {
//        file << ::std::setprecision(6) << volSurf.maturities[i];
//
//        for (size_t j = 0; j < volSurf.moneyness.size(); ++j)
//        {
//            file << "," << ::std::setprecision(6) << volSurf.vols[i][j];
//        }
//
//        file << "\n";
//    }
//
//    file.close();
//}