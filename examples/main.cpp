#include "Utils/Data/CSV/CSVRead.hpp"
#include "Utils/Log.hpp"
#include "Utils/StopWatch/StopWatch.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Models/LocalVol/LocalVol.hpp"
#include "Models/Heston/HestonPricer/HestonPricer.hpp"
#include "Models/Heston/HestonConfig.hpp"
#include "Models/Heston/HestonCalibrator/HestonCalibrator.hpp"
#include "Errors/Errors.hpp"  
#include "Math/MathFunctions/MathFunctions.hpp"
#include "Math/Quadrature/TanHSinH.hpp"
#include "Math/Calibration/Ceres/CeresPolicy.hpp"
#include "Math/Calibration/Ceres/CeresConfig.hpp"
#include "Math/Interpolation/Interpolation.hpp"


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
#include <limits>
#include <memory>

using namespace uv;

int main(int argc, char* argv[])
{
    try
    {
        // Choose the CSV path:
        //  - If the user passed a command-line argument, use that path (argv[1]).
        //  - Otherwise, fall back to your default CSV inside the repo.
        const std::filesystem::path path = (argc > 1) ?
            std::filesystem::path{ argv[1] } :
            std::filesystem::path{ "data/inputs/VolSurface_SPY_04072011.csv"};

        // Set logger
        UV_LOG_TO_FILE("calibration.log");
        UV_LOG_CONSOLE(true);

        // Start timer
        StopWatch timer;
        timer.StartStopWatch();

        // Define market data
        MarketData mktData
        { 0.0,          // r
          0.0,          // q
          485.77548     // S
        };

        // Generate volatility surface instance
        VolSurface mktVolSurf{ readVolSurface(path.string(), mktData) };
        //mktVolSurf.printTotVar();

        // Initialize NLopt Calibrator instance
        CalibratorNLopt<5, nlopt::LD_SLSQP> nloptOptimizer
        {
            NLoptConfig<5>
            {
                1e-12,                             // eps
                1e-9,                              // tol
                1e-10,                             // ftolRel
                10000,                             // maxEval
                { "a", "b", "rho", "m", "sigma" }  // Parameter names 
            }
        };

        // Calibrate SVI surface
        auto [sviSlices, sviVolSurface] = svi::calibrate(mktVolSurf, nloptOptimizer);
        sviVolSurface.printVol();

        // Build the Local Volatility surface using SVI parameters
        const VolSurface lvVolSurface{local_vol::build(sviVolSurface, sviSlices)};
        lvVolSurface.printVol();

        // Initialize integration quadrature
        static constexpr std::size_t HestonNodes = 300;
        const TanHSinH<HestonNodes> quad{};

        // Initialize HestonPricer instance
        HestonPricer hestonPricer
        {
            std::make_shared<const TanHSinH<HestonNodes>>(quad),
            {
                -2.0,
                2.0
            }
        };

        // Initialize Ceres calibrator instance
        CalibratorCeres
            <
                5,                                // Number of calibration parameters (kappa, theta, sigma, rho, v0)
                CeresPolicy
                  <
                    ceres::HuberLoss,            // Robust loss function type (Huber loss for outlier resistance)
                    ceres::LEVENBERG_MARQUARDT,  // Trust region strategy (LM algorithm)
                    ceres::DENSE_QR              // Linear solver type (dense QR decomposition for small problems)
                >
            > ceresOptimizer
        { 
            CeresConfig<5>
            {   
                
                1000,                                            // Max evaluations
                1e-16,                                           // Function tolerance
                1e-16,                                           // Parameter tolerance
                { "kappa", "theta", "sigma", "rho", "v0" },      // Parameter names
                1e-16,                                           // Gradient tolerance
                1.0,                                             // Loss scaling parameter   
                false                                            // Logs the full Ceres calibration report 
            }
        };

        // Calibrate the Heston model
        VolSurface hestonVolurface
        { 
            heston_calibrator::calibrate
            (
                sviVolSurface,
                hestonPricer,
                ceresOptimizer
            )
        };
        hestonVolurface.printBSCall();

        // End and log timer
        timer.LogTime<std::milli>();

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








//{
//    const double kappa = 5.0;
//    const double theta = 0.045324;
//    const double sigma = 2.606662;
//    const double rho = -0.759016;
//    const double v0 = 0.371132;
//    double S = 485.77548;
//    double T = 3.0;
//    const double r = 0.0;
//    const double q = 0.0;
//
//    const SliceData& sliceData = sviVolSurface.slices().back();
//    const std::vector<double>& strikes = sliceData.K();
//    //std::cout << std::fixed << std::setprecision(17);
//
//    for (std::size_t i = 0; i < strikes.size(); ++i)
//    {
//        std::cout << "K: " << sliceData.mny()[i] << " ";
//        std::cout << hestonPricer.callPrice(
//            kappa,
//            theta,
//            sigma,
//            rho,
//            v0,
//            T,
//            S * std::exp((r - q) * T),
//            r,
//            strikes[i]) << "\n";
//    }
//
//        }
//
//
//// Parameters where Heston = BS
//const double kappa = 1.0;     // irrelevant if sigma = 0
//const double theta = 0.04;    // variance = 0.04 -> vol = 0.2
//const double sigma = 1e-15;    // almost zero -> constant variance
//
//const double rho = 0.0;
//const double v0 = theta;   // constant variance
//double S = 100.0;
//const double K = 100.0;
//double T = 1.0;
//const double r = 0.02;
//const double q = 0.0;
//
//double hestonPrice;
//
//for (int i = 0; i < 2; ++i)
//{
//    hestonPrice = hestonPricer.callPrice(
//        kappa,
//        theta,
//        sigma,
//        rho,
//        v0,
//        T,
//        S * std::exp((r - q) * T),
//        r,
//        K
//    );
//}
//
//
//
//const double bsPrice = blackScholes(T, r, q, std::sqrt(theta), S, K);
//std::cout << std::setprecision(17);
//std::cout << "Heston price: " << hestonPrice << "\n";
//std::cout << "BS price:     " << bsPrice << "\n";
//std::cout << "Abs diff:     " << std::fabs(hestonPrice - bsPrice) << "\n";
//
//// simple check
//if (std::fabs(hestonPrice - bsPrice) < 1e-5)
//    std::cout << "✅ Heston matches Black–Scholes.\n";
//else
//std::cout << "❌ Mismatch — check integrand or scaling.\n";
//
//
//T = 0.0;
//S = 105.0;
//
//
//hestonPrice = hestonPricer.callPrice(
//    kappa,
//    theta,
//    sigma,
//    rho,
//    v0,
//    T,
//    S * std::exp((r - q) * T),
//    r,
//    K
//);
//
//const double intrinsic = std::max(S - K, 0.0);
//const double diff = std::abs(hestonPrice - intrinsic);
//
//std::cout << "\n--- T → 0 Limit Test ---\n";
//std::cout << "Heston Price  : " << hestonPrice << '\n';
//std::cout << "Intrinsic Val.: " << intrinsic << '\n';
//std::cout << "Abs Difference: " << diff << '\n';
//
//if (diff < 1e-6)
//    std::cout << "✅ Passed: Heston price converges to intrinsic value as T → 0.\n";
//else
//std::cout << "❌ Failed: Heston price deviates from intrinsic value.\n";
