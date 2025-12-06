#include "Utils/IO/CSVRead.hpp"
#include "Utils/IO/Log.hpp"
#include "Utils/Aux/Helpers.hpp"
#include "Utils/Aux/StopWatch.hpp"
#include "Core/VolSurface.hpp"
#include "Core/MarketData.hpp"
#include "Models/SVI/SVI.hpp"
#include "Models/LocalVol/LocalVol.hpp"
#include "Models/Heston/Pricer.hpp"
#include "Models/Heston/Config.hpp"
#include "Models/Heston/Calibrator.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/Types.hpp"
#include "Math/Functions.hpp"
#include "Math/Quadratures/TanHSinH.hpp"
#include "Calibration/Ceres/Policy.hpp"
#include "Calibration/Ceres/Config.hpp"
#include "Math/Interpolation.hpp"


#include <chrono>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <string_view>
#include <numeric>
#include <cassert>
#include <utility>
#include <numbers>
#include <limits>
#include <memory>


using namespace uv;
using namespace models;
using namespace utils;
using namespace core;
using namespace math;


int main(int argc, char* argv[])
{
    try
    {
        // ---------- Configurations ----------
        
        // Choose the CSV path
        const std::filesystem::path path = (argc > 1) ?
            std::filesystem::path{ argv[1] } :
            std::filesystem::path{ "data/inputs/VolSurface_SPY_04072011.csv"};

        // Set logger
        UV_LOG_TO_FILE("calibration.log");
        UV_LOG_CONSOLE(true);

        // Start timer
        StopWatch timer;
        timer.StartStopWatch();

        // ---------- Market data ----------

        MarketData mktData
        { 
          Real(0.0),          // r
          Real(0.0),          // q
          Real(485.77548)     // S
        };

        VolSurface mktVolSurf{ readVolSurface(path.string(), mktData) };
        //mktVolSurf.printTotVar(); 

        // ---------- SVI Calibration ----------

        cal::nlopt::Calibrator<5, nlopt::LD_SLSQP> nloptCalibrator
        {
            cal::nlopt::Config<5>
            {
                1e-12,                            
                1e-9,                              
                1e-10,                            
                10000,                             
                { "a", "b", "rho", "m", "sigma" }   
            }
        };

        auto [sviSlices, sviVolSurface] = svi::calibrate(mktVolSurf, nloptCalibrator);
        sviVolSurface.printVol();

        // ---------- Build Local Volatility surface ----------

        const VolSurface lvVolSurface{localvol::buildSurface(sviVolSurface, sviSlices)};
        lvVolSurface.printVol();

        // BP()

        // ---------- Heston model calibration ----------

        static constexpr std::size_t HestonNodes = 300;
        const TanHSinH<HestonNodes> quad{};
        
        heston::Pricer hestonPricer
        {
            std::make_shared<const TanHSinH<HestonNodes>>(quad),
            {
                Real(-2.0),   // Damping parameter ITM
                Real(2.0)     // Damping parameter OTM
            }
        };

        cal::ceres::Calibrator<5, cal::ceres::Policy
                <
                    ceres::HuberLoss,             
                    ceres::LEVENBERG_MARQUARDT,   
                    ceres::DENSE_QR               
                >
            > ceresOptimizer
        { 
            cal::ceres::Config<5>
            {   
                
                1000,                                          
                1e-16,                                        
                1e-16,                                         
                { "kappa", "theta", "sigma", "rho", "v0" },      
                1e-16,                                           
                1.0,                                               
                false                                           
            }
        };

        VolSurface hestonVolurface
        { 
            heston::calibrator::calibrate
            (
                sviVolSurface,
                hestonPricer,
                ceresOptimizer
            )
        };

        hestonVolurface.printVol();

        // ---------- Outputs ----------

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
//    const Real kappa = Real(5.0);
//    const Real theta = 0.045324;
//    const Real sigma = 2.606662;
//    const Real rho = -0.759016;
//    const Real v0 = 0.371132;
//    Real S = 485.77548;
//    Real T = Real(3.0);
//    const Real r = Real(0.0);
//    const Real q = Real(0.0);
//
//    const SliceData& sliceData = sviVolSurface.slices().back();
//    const Vector<Real>& strikes = sliceData.K();
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
//const Real kappa = Real(1.0);     // irrelevant if sigma = 0
//const Real theta = 0.04;    // variance = 0.04 -> vol = 0.2
//const Real sigma = 1e-15;    // almost zero -> constant variance
//
//const Real rho = Real(0.0);
//const Real v0 = theta;   // constant variance
//Real S = Real(100.0);
//const Real K = Real(100.0);
//Real T = Real(1.0);
//const Real r = 0.02;
//const Real q = Real(0.0);
//
//Real hestonPrice;
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
//const Real bsPrice = blackScholes(T, r, q, std::sqrt(theta), S, K);
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
//T = Real(0.0);
//S = Real(105.0);
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
//const Real intrinsic = std::max(S - K, Real(0.0));
//const Real diff = std::abs(hestonPrice - intrinsic);
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
