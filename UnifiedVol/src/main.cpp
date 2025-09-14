#include "Utils/CSVRead.hpp"
#include "Utils/CSVWrite.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"

#include <iostream>

int main() 
{
    MarketData mkt{0.0, 0.0};

    VolSurface mktVolSurf{ readVolSurface("data/inputs/VolSurface_SPY_04072011.csv", mkt) };


    mktVolSurf.printImpVol();
    std::cout << "\n";
    mktVolSurf.printTotVar();
    std::cout << "\n";

    SVI svi{ mktVolSurf };

    //svi.modelVolSurf.printConsole();

    //writeVolSurface(svi.modelVolSurf, "data/SVI_VolSurface_SPY_04072011.csv");

    return 0;
 }