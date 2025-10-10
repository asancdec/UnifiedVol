#include "Utils/Data/CSVRead.hpp"
#include "Utils/Log.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Errors/Errors.hpp"  
#include "Quadrature/StandardGL.hpp"

#include <chrono>
#include <filesystem>
#include <format>
#include <iostream>
#include <iomanip>
#include <string_view>
#include <vector>


int main(int argc, char* argv[])
{
    try 
    {   

        // Choose the CSV path:
        //  - If the user passed a command-line argument, use that path (argv[1]).
        //  - Otherwise, fall back to your default CSV inside the repo.
        const std::filesystem::path path = (argc > 1) ? 
            std::filesystem::path{ argv[1] }:   
            std::filesystem::path{ "data/inputs/VolSurface_SPY_04072011.csv" };

        UV_LOG_TO_FILE("calibration.log");
        UV_LOG_CONSOLE(true);

        // Start timer
        const auto t0{ std::chrono::high_resolution_clock::now() };

        MarketData mkt{ 0.0, 0.0 };

        VolSurface mktVolSurf{ readVolSurface(path.string(), mkt) };
        mktVolSurf.printTotVar();

        bool isValidateResults{ true };
        SVI svi{ mktVolSurf, isValidateResults };

        VolSurface sviVolSurf{ svi.getVolSurf() };

        sviVolSurf.printVol();

        StandardGL stgl{ 64 };
        stgl.printGrid();

        // End timer
        const auto t1 = std::chrono::high_resolution_clock::now();
        const auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
        UV_INFO(std::format("Done in {} ms", dt.count()));

        return EXIT_SUCCESS;
    }
    catch (const uv::UnifiedVolError& e)
    {
        std::cerr << e.what() << '\n';
        return 2; // domain error
    }
    catch (const std::exception& e)
    {
        std::cerr << "std::exception: " << e.what() << '\n';
        return 1; // generic error
    }
}
