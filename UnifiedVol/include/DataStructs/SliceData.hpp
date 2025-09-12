/**
* RawSVIParams.hpp
* Author: Alvaro Sanchez de Carlos
* Date: 08/18/2025
*/

#ifndef SLICE_DATA_HPP
#define SLICE_DATA_HPP

#include <vector>

// Struct to hold market data for one maturity
struct SliceData 
{
    std::vector<double> k;    // log-forward moneyness
    std::vector<double> wMkt; // total implied variance (sigma^2 * T)
};

#endif // SLICE_DATA_HPP
