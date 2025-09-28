#include "Utils/Data/CSVRead.hpp"
#include "Utils/Data/CSVWrite.hpp"
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

        // Start timer
        const auto t0{ std::chrono::high_resolution_clock::now() };

        //MarketData mkt{ 0.0, 0.0 };

        //VolSurface mktVolSurf{ readVolSurface(path.string(), mkt) };
        //mktVolSurf.printVol();
        //std::cout << "\n";

        //bool isPrintResults{ true };
        //SVI svi{ mktVolSurf, isPrintResults };
        //std::cout << "\n";

        //VolSurface sviVolSurf{ svi.getVolSurf() };

        //sviVolSurf.printTotVar();

        StandardGL stgl{ 64, 0.0 };
        stgl.printGrid();

        // 1) f(x) = 1  → ∫_0^∞ e^{-x} dx = 1
        std::cout << "I0  = " << stgl.eval([](double) { return 1.0; })
            << "   (expected 1)\n";

        // 2) f(x) = x  → ∫_0^∞ x e^{-x} dx = 1!
        std::cout << "I1  = " << stgl.eval([](double x) { return x; })
            << "   (expected 1)\n";

        // 3) f(x) = x^2 → ∫_0^∞ x^2 e^{-x} dx = 2!
        std::cout << "I2  = " << stgl.eval([](double x) { return x * x; })
            << "   (expected 2)\n";

        // 4) f(x) = x^3 → ∫_0^∞ x^3 e^{-x} dx = 3! = 6
        std::cout << "I3  = " << stgl.eval([](double x) { return x * x * x; })
            << "   (expected 6)\n";

        // 5) f(x) = e^{-0.5x} → ∫ e^{-x} e^{-0.5x} dx = ∫ e^{-1.5x} dx = 1/1.5 = 2/3
        std::cout << "Ie  = " << stgl.eval([](double x) { return std::exp(-0.5 * x); })
            << "   (expected ~0.6666666666667)\n";

        // 6) f(x) = cos(x) → ∫ e^{-x} cos(x) dx = a/(a^2+b^2) with a=1,b=1 → 1/2
        std::cout << "Icos= " << stgl.eval([](double x) { return std::cos(x); })
            << "   (expected 0.5)\n";

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
