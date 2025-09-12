#include "Utils/CSVRead.hpp"
#include "Utils/CSVWrite.hpp"
#include "DataStructs/VolSurface.hpp"
#include "VolModels/SVI.hpp"


int main() 
{

    VolSurface impVolSurf{ readVolSurface("data/VolSurface_SPY_04072011.csv") };

    impVolSurf.printConsole();

    //SVI svi{impVolSurf};

    //svi.modelVolSurf.printConsole();

    //writeVolSurface(svi.modelVolSurf, "data/SVI_VolSurface_SPY_04072011.csv");

    return 0;
 }