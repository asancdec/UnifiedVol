#include "Utils/Data/CSVRead.hpp"
#include "Utils/Data/CSVWrite.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Errors/Errors.hpp"  

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

        // Start timer
        const auto t0{ std::chrono::high_resolution_clock::now() };

        MarketData mkt{ 0.0, 0.0 };

        VolSurface mktVolSurf{ readVolSurface(path.string(), mkt) };
        mktVolSurf.printVol();
        std::cout << "\n";

        bool isPrintResults{ true };
        SVI svi{ mktVolSurf, isPrintResults };
        std::cout << "\n";

        VolSurface sviVolSurf{ svi.getVolSurf() };

        sviVolSurf.printVol();
        std::cout << "\n";

        // End timer
        const auto t1{ std::chrono::high_resolution_clock::now() };
        const auto dt{ std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0) };
        std::cout << std::format("Done in {} ms\n", dt.count());

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
