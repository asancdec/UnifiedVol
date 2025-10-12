/**
* Heston.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/Heston/Heston.hpp"
#include <array>

Heston::Heston(const VolSurface& volSurf) : calVolSurf_(volSurf)
{

}



//std::array<double, 5> Heston::initGuess() noexcept
//{
//    const double b{ 0.1 };
//    const double rho{ -0.5 };
//    const double m{ 0.1 };
//    const double sigma{ 0.1 };
//    const double a{ slice.atmWT() - b * (-rho * m + std::hypot(m, sigma)) };
//
//    return { a, b, rho, m, sigma };
//}
//
//std::array<double, 5> Heston::lowerBounds() noexcept
//{
//    return
//    {
//        0.001,                     // kappa
//        0.001,                     // theta            
//        0.001,                     // sigma
//        -0.999                     // rho
//        0.001                      // v0
//    };
//}
//
//std::array<double, 5> Heston::upperBounds() noexcept
//{
//    return
//    {
//        0.001,                     // kappa
//        0.001,                     // theta            
//        0.001,                     // sigma
//        0.999                      // rho
//        0.001                      // v0
//    };
//}

