#include "Utils/Data/CSVRead.hpp"
#include "Utils/Log.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Models/Heston/HestonPricer.hpp"
#include "Models/Heston/HestonConfig.hpp"
#include "Models/Heston/HestonCalibrator.hpp"
#include "Errors/Errors.hpp"  
#include "Math/MathFunctions/MathFunctions.hpp"
#include "Math/Quadrature/TanHSinH.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <string_view>
#include <vector>
#include <numeric>
#include <cassert>
#include <utility>
#include <numbers>
#include <memory>

using namespace uv;

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

        MarketData mktData
        { 0.0,          // r
          0.0,          // q
          485.77548     // S
        };

        VolSurface mktVolSurf{ readVolSurface(path.string(), mktData) };
        mktVolSurf.printVol();

        const auto t0{ std::chrono::high_resolution_clock::now() };

        for (auto& slice : mktVolSurf.slices())
        {
            auto copy = slice.callBS();
            slice.setCallBS(copy);
        }


        const auto t1 = std::chrono::high_resolution_clock::now();
        mktVolSurf.printVol();


        CalibratorConfig<5> sviConfig
        {
            1e-12,                             // eps
            1e-9,                              // tol
            1e-10,                             // ftolRel
            10000,                             // maxEval
            { "a", "b", "rho", "m", "sigma" }  // Parameter names
        };

        // Initialize Calibrator instance
        CalibratorNLopt<5> nloptOptimizer
        {
            sviConfig,
            nlopt::LD_SLSQP
        };

        SVI svi{};

        SVIReport sviReport{ svi.calibrate(mktVolSurf, nloptOptimizer) };

        TanHSinH quad{ 1000 };


        HestonPricer hestonPricer
        {
            std::make_shared<const TanHSinH>(quad),
            {
                -3.5,
                3.5
            }
        };



        // Parameters where Heston = BS
        const double kappa = 1.0;     // irrelevant if sigma = 0
        const double theta = 0.04;    // variance = 0.04 -> vol = 0.2
        const double sigma = 1e-15;    // almost zero -> constant variance

        const double rho = 0.0;
        const double v0 = theta;   // constant variance
        double S = 100.0;
        const double K = 100.0;
        double T = 1.0;
        const double r = 0.02;
        const double q = 0.0;


        double hestonPrice = hestonPricer.callPrice(
            kappa,
            theta,
            sigma,
            rho,
            v0,
            T,
            S * std::exp((r - q) * T),
            r,
            K
        );


        const double bsPrice = blackScholesCall(S, K, r, q, T, std::sqrt(theta));
        std::cout << std::setprecision(17);
        std::cout << "Heston price: " << hestonPrice << "\n";
        std::cout << "BS price:     " << bsPrice << "\n";
        std::cout << "Abs diff:     " << std::fabs(hestonPrice - bsPrice) << "\n";

        // simple check
        if (std::fabs(hestonPrice - bsPrice) < 1e-5)
            std::cout << "✅ Heston matches Black–Scholes.\n";
        else
            std::cout << "❌ Mismatch — check integrand or scaling.\n";


        T = 0.0;
        S = 105.0;


        hestonPrice = hestonPricer.callPrice(
            kappa,
            theta,
            sigma,
            rho,
            v0,
            T,
            S * std::exp((r - q) * T),
            r,
            K
        );

        const double intrinsic = std::max(S - K, 0.0);
        const double diff = std::abs(hestonPrice - intrinsic);

        std::cout << "\n--- T → 0 Limit Test ---\n";
        std::cout << "Heston Price  : " << hestonPrice << '\n';
        std::cout << "Intrinsic Val.: " << intrinsic << '\n';
        std::cout << "Abs Difference: " << diff << '\n';

        if (diff < 1e-6)
            std::cout << "✅ Passed: Heston price converges to intrinsic value as T → 0.\n";
        else
            std::cout << "❌ Failed: Heston price deviates from intrinsic value.\n";


        CalibratorConfig<5> hestonConfig
        {
            1e-12,                                      // eps
            1e-9,                                       // tol
            1e-10,                                      // ftolRel
            10000,                                      // maxEval
            { "kappa", "theta", "sigma", "rho", "v0" }  // Parameter names
        };

        CalibratorCeres<5> ceresOptimizer{ hestonConfig };


        HestonCalibrator hestonCalibrator{};
        //hestonCalibrator.calibrate
        //(
        //    mktVolSurf,
        //    std::make_shared<const HestonPricer>(hestonPricer),
        //    ceresOptimizer
        //);


        // End time

        const auto dt = std::chrono::duration_cast<::std::chrono::milliseconds>(t1 - t0);
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


