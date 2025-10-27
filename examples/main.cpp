#include "Utils/Data/CSVRead.hpp"
#include "Utils/Log.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Models/Heston/Heston.hpp"
#include "Errors/Errors.hpp"  
#include "Math/Quadrature/TanHSinH.hpp"
#include <chrono>
#include <filesystem>
#include <format>
#include <iostream>
#include <iomanip>
#include <string_view>
#include <vector>
#include <numeric>
#include <cassert>
#include <utility>
#include <numbers>


inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double blackScholesCall(double S, double K, double r, double q, double T, double vol) {
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * vol * vol) * T) / (vol * sqrtT);
    const double d2 = d1 - vol * sqrtT;
    return S * std::exp(-q * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
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
        
        TanHSinH quad{ 1000 };
        quad.printGrid();

        std::cout << std::setprecision(17);


        const double PI = std::numbers::pi;
        const double EULERG = std::numbers::egamma;; // Euler–Mascheroni γ

        // 1) ∫_0^∞ e^{-z} dz = 1
        auto f1 = [](double z) { return std::exp(-z); };
        double v1 = quad.integrateZeroToInf(f1);
        std::cout << "I1 = ∫_0^∞ e^{-z} dz            : value = " << v1
            << "   exact = " << 1.0
            << "   abs err = " << std::fabs(v1 - 1.0) << std::endl;

        // 2) ∫_0^∞ 1/(1+z^2) dz = π/2
        auto f2 = [](double z) { return 1.0 / (1.0 + z * z); };
        double v2 = quad.integrateZeroToInf(f2);
        std::cout << "I2 = ∫_0^∞ 1/(1+z^2) dz        : value = " << v2
            << "   exact = " << (PI / 2.0)
            << "   abs err = " << std::fabs(v2 - (PI / 2.0)) << std::endl;

        // 3) ∫_0^∞ z^2 e^{-z} dz = Γ(3) = 2
        auto f3 = [](double z) { return z * z * std::exp(-z); };
        double v3 = quad.integrateZeroToInf(f3);
        std::cout << "I3 = ∫_0^∞ z^2 e^{-z} dz       : value = " << v3
            << "   exact = " << 2.0
            << "   abs err = " << std::fabs(v3 - 2.0) << std::endl;

        // 4) ∫_0^∞ e^{-z} log(z) dz = -γ
        auto f4 = [](double z) { return std::exp(-z) * std::log(z); };
        double v4 = quad.integrateZeroToInf(f4);
        std::cout << "I4 = ∫_0^∞ e^{-z} log(z) dz    : value = " << v4
            << "   exact = " << (-EULERG)
            << "   abs err = " << std::fabs(v4 + EULERG) << std::endl;

        // 5) ∫_0^∞ √z e^{-z} dz = Γ(3/2) = √π / 2
        auto f5 = [](double z) { return std::sqrt(z) * std::exp(-z); };
        double v5 = quad.integrateZeroToInf(f5);
        std::cout << "I5 = ∫_0^∞ √z e^{-z} dz        : value = " << v5
            << "   exact = " << (std::sqrt(PI) / 2.0)
            << "   abs err = " << std::fabs(v5 - (std::sqrt(PI) / 2.0)) << std::endl;

        // 6) ∫_0^∞ e^{-2z} dz = 1/2
        auto f6 = [](double z) { return std::exp(-2.0 * z); };
        double v6 = quad.integrateZeroToInf(f6);
        std::cout << "I6 = ∫_0^∞ e^{-2z} dz          : value = " << v6
            << "   exact = " << 0.5
            << "   abs err = " << std::fabs(v6 - 0.5) << std::endl;

        // 7) ∫_0^∞ z e^{-3z} dz = 1/9
        auto f7 = [](double z) { return z * std::exp(-3.0 * z); };
        double v7 = quad.integrateZeroToInf(f7);
        std::cout << "I7 = ∫_0^∞ z e^{-3z} dz        : value = " << v7
            << "   exact = " << (1.0 / 9.0)
            << "   abs err = " << std::fabs(v7 - (1.0 / 9.0)) << std::endl;

        // 8) ∫_0^∞ z^4 e^{-z} dz = Γ(5) = 4! = 24
        auto f8 = [](double z) { return (z * z * z * z) * std::exp(-z); };
        double v8 = quad.integrateZeroToInf(f8);
        std::cout << "I8 = ∫_0^∞ z^4 e^{-z} dz       : value = " << v8
            << "   exact = " << 24.0
            << "   abs err = " << std::fabs(v8 - 24.0) << std::endl;

        // 9) ∫_0^∞ 1/(1+z^3) dz = 2π / (3√3)
        auto f9 = [](double z) { return 1.0 / (1.0 + z * z * z); };
        double v9 = quad.integrateZeroToInf(f9);
        const double exact9 = 2.0 * PI / (3.0 * std::sqrt(3.0));
        std::cout << "I9 = ∫_0^∞ 1/(1+z^3) dz        : value = " << v9
            << "   exact = " << exact9
            << "   abs err = " << std::fabs(v9 - exact9) << std::endl;

        // 10) ∫_0^∞ 1/(1+z^4) dz = π / (2√2)
        auto f10 = [](double z) { return 1.0 / (1.0 + (z * z * z * z)); };
        double v10 = quad.integrateZeroToInf(f10);
        const double exact10 = PI / (2.0 * std::sqrt(2.0));
        std::cout << "I10 = ∫_0^∞ 1/(1+z^4) dz       : value = " << v10
            << "   exact = " << exact10
            << "   abs err = " << std::fabs(v10 - exact10) << std::endl;

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




//Heston heston(sviVolSurf, gl);

//// Parameters where Heston = BS
//const double kappa = 1.0;     // irrelevant if sigma = 0
//const double theta = 0.04;    // variance = 0.04 -> vol = 0.2
//const double sigma = 1e-8;    // almost zero -> constant variance
//const double rho = 0.0;
//const double v0 = theta;   // constant variance
//const double S = 100.0;
//const double K = 100.0;
//const double T = 1.0;
//const double r = 0.02;
//const double q = 0.0;

//const double hestonPrice = heston.callPrice(kappa, theta, sigma, rho, v0, T, S, r, q, K);
//const double bsPrice = blackScholesCall(S, K, r, q, T, std::sqrt(theta));

//std::cout << "Heston price: " << hestonPrice << "\n";
//std::cout << "BS price:     " << bsPrice << "\n";
//std::cout << "Abs diff:     " << std::fabs(hestonPrice - bsPrice) << "\n";

//// simple check
//if (std::fabs(hestonPrice - bsPrice) < 1e-5)
//    std::cout << "✅ Heston matches Black–Scholes.\n";
//else
//    std::cout << "❌ Mismatch — check integrand or scaling.\n";