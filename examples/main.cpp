#include "Utils/Data/CSVRead.hpp"
#include "Utils/Log.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Errors/Errors.hpp"  
#include "Math/Quadrature/GaussLaguerre.hpp"

#include <chrono>
#include <filesystem>
#include <format>
#include <iostream>
#include <iomanip>
#include <string_view>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>


static inline double relerr(double a, double b) {
    if (b == 0.0) return (a == 0.0) ? 0.0 : std::numeric_limits<double>::infinity();
    return std::abs((a - b) / b);
}

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

        MarketData mkt
        { 0.0,          // r
          0.0,          // q
          485.77548     // S
        };

        VolSurface mktVolSurf{ readVolSurface(path.string(), mkt) };

        SVI svi{ mktVolSurf  };

        VolSurface sviVolSurf{ svi.getVolSurf() };
        //sviVolSurf.printBSCall();
        

        
        GaussLaguerre gl{ 64, 0.0};
        gl.printGrid();

        // Plain integrals (unweighted):
        // f0(x) = e^{-1.5 x}            ⇒ ∫_0^∞ f0 = 1/1.5
        // f1(x) = e^{-2.0 x} cos(0.8 x) ⇒ ∫_0^∞ f1 = a/(a^2 + b^2) with a=2, b=0.8
        // f2(x) = x e^{-2.0 x}          ⇒ ∫_0^∞ f2 = 1/a^2 with a=2
        auto g0 = [&](double x) { return std::exp(-1.5 * x); };
        auto g1 = [&](double x) { return std::exp(-2.0 * x) * std::cos(0.8 * x); };
        auto g2 = [&](double x) { return x * std::exp(-2.0 * x); };

        const double E0 = 1.0 / 1.5;                       // 0.666666666666...
        const double E1 = 2.0 / (2.0 * 2.0 + 0.8 * 0.8);       // 2 / (4 + 0.64) = 0.431034482758...
        const double E2 = 1.0 / (2.0 * 2.0);                 // 0.25

        const double I0 = gl.evalUnweighted(g0);
        const double I1 = gl.evalUnweighted(g1);
        const double I2 = gl.evalUnweighted(g2);

        auto show = [](const char* name, double I, double E) {
            std::cout << name << "  I=" << I
                << "  exact=" << E
                << "  rel=" << relerr(I, E) << "\n";
            };

        std::cout << "\n[Unweighted integrals, α=0]\n";
        show("f0: e^{-1.5 x}               ", I0, E0);
        show("f1: e^{-2 x} cos(0.8 x)      ", I1, E1);
        show("f2: x e^{-2 x}               ", I2, E2);


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
