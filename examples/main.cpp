// SPDX-License-Identifier: Apache-2.0
/*
 * File:        main.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#include "Utils/PCH.hpp"


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
            std::filesystem::path{ "data/VolSurface_SPY_04072011.csv"};

        // Set logger
        UV_LOG_TO_FILE("calibration.log");
        UV_LOG_CONSOLE(true);

        // Start timer
        StopWatch timer;
        timer.StartStopWatch();

        // ---------- Market data ----------

        MarketData marketData
        { 
          Real(0.0),          // r
          Real(0.0),          // q
          Real(485.77548)     // S
        };

        VolSurface mktVolSurface{ readVolSurface(path.string(), marketData) };
        mktVolSurface.printTotVar();

        // ---------- SVI Calibration ----------

        opt::nlopt::Optimizer<4, nlopt::LD_SLSQP> nloptOptimizer
        {
            opt::nlopt::Config<4>
            {
                1e-12,                            
                1e-9,                              
                1e-10,                            
                10000,                             
                { "b", "rho", "m", "sigma" }   
            }
        };

        Vector<svi::Params> sviParams
        { 
            svi::calibrate
            (
                mktVolSurface.tenors(),
                mktVolSurface.logKFMatrix(),
                mktVolSurface.totVarMatrix(),
                nloptOptimizer
            ) 
        
        };

        const VolSurface sviVolSurface
        {
            svi::buildSurface
            (
                mktVolSurface,
                sviParams
            )
        };
       
        sviVolSurface.printTotVar();

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

        opt::ceres::Optimizer<5, opt::ceres::Policy
                <
                    ceres::HuberLoss,             
                    ceres::LEVENBERG_MARQUARDT,   
                    ceres::DENSE_QR               
                >
            > ceresOptimizer
        { 
            opt::ceres::Config<5>
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

        hestonPricer.setParams
        (
            heston::calibrator::calibrate
            (
                sviVolSurface.tenors(),
                sviVolSurface.strikes(),
                sviVolSurface.forwards(),
                sviVolSurface.rates(),
                sviVolSurface.callPrices(),
                hestonPricer,
                ceresOptimizer
            )
        );

        const VolSurface hestonVolSurface
        { 
            heston::calibrator::buildSurface
            (
                sviVolSurface, 
                hestonPricer
            ) 
        };

        hestonVolSurface.printVol();

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