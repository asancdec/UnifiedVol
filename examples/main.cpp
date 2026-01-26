// SPDX-License-Identifier: Apache-2.0
/*
 * File:        main.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Example script of model calibrations and pricing.
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
        // ---------- Configure ----------

        // Choose the CSV path
        const std::filesystem::path path =
            (argc > 1) ? std::filesystem::path{argv[1]}
                       : std::filesystem::path{"data/VolSurface_SPY_04072011.csv"};

        // Set logger
        UV_LOG_TO_FILE("calibration.log");
        UV_LOG_CONSOLE(true);

        // Start timer
        StopWatch timer;
        timer.StartStopWatch();

        // ---------- Read ----------

        MarketData<Real> marketData{.r = 0.00, .q = 0.00, .S = 485.77548};

        VolSurface<Real> mktVolSurface{readVolSurface<Real>(path.string(), marketData)};

        mktVolSurface.printTotVar();

        // ---------- SVI calibration ----------

        opt::nlopt::Optimizer<4, nlopt::LD_SLSQP> nloptOptimizer{opt::nlopt::Config<4>{
            .eps = 1e-12,
            .tol = 1e-12,
            .ftolRel = 1e-12,
            .maxEval = 10000,
            .paramNames = {"b", "rho", "m", "sigma"}}};

        Vector<svi::Params<Real>> sviParams{svi::calibrate<Real>(
            mktVolSurface.tenors(),
            mktVolSurface.logKFMatrix(),
            mktVolSurface.totVarMatrix(),
            nloptOptimizer
        )};

        const VolSurface<Real> sviVolSurface{
            svi::buildSurface<Real>(mktVolSurface, sviParams)};

        sviVolSurface.printBSCall();

        // ---------- Local volatility calibration ----------

        constexpr std::size_t nT{100}; // Time steps
        constexpr std::size_t nX{20};  // Space steps
        constexpr Real xBound{2.5};    // log(K/F) grid bound

        // Declare payoff function
        // TODO: Parametrize as a function
        auto payoff = [](Real S, Real K) -> Real
        {
            return (S > K) ? (S - K) : 0.0;
        };

        // Allocate to heap
        auto lvPricer = std::make_unique<
            localvol::Pricer<Real, nT, nX, math::interp::PchipInterpolator<Real>>>(
            payoff,
            xBound
        );

        localvol::Surface<Real> localVolSurface{
            localvol::calibrator::
                calibrate<Real, nT, nX, math::interp::PchipInterpolator<Real>>(
                    sviVolSurface.callPrices(),
                    sviVolSurface.tenors(),
                    sviVolSurface.logKFMatrix(),
                    sviVolSurface.totVarMatrix(),
                    *lvPricer
                )};

        // ---------- Heston model calibration ----------

        constexpr std::size_t HestonNodes{300};
        const integration::TanHSinH<Real, HestonNodes> quad{};

        heston::Pricer<Real, HestonNodes> hestonPricer{
            std::make_shared<const integration::TanHSinH<Real, HestonNodes>>(quad),
            {
                -2.0, // Damping parameter ITM
                2.0   // Damping parameter OTM
            }};

        opt::ceres::Optimizer<
            5,
            opt::ceres::Policy<
                void, // No Robust Loss
                ceres::LEVENBERG_MARQUARDT,
                ceres::DENSE_NORMAL_CHOLESKY>>
            lmOptimizer{opt::ceres::Config<5>{

                .maxEval = 1000,
                .functionTol = 1e-12,
                .paramTol = 1e-12,
                .gradientTol = 1e-12,
                .paramNames = {"kappa", "theta", "sigma", "rho", "v0"},
                .verbose = false}};

        hestonPricer.setParams(heston::calibrator::calibrate<Real>(
            sviVolSurface.tenors(),
            sviVolSurface.strikes(),
            sviVolSurface.forwards(),
            sviVolSurface.rates(),
            sviVolSurface.callPrices(),
            hestonPricer,
            lmOptimizer,
            opt::WeightATM{.wATM = 8.0, .k0 = 0.3}
        ));

        const VolSurface<Real> hestonVolSurface{
            heston::calibrator::buildSurface<Real>(sviVolSurface, hestonPricer)};

        hestonVolSurface.printBSCall();

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