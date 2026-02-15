// SPDX-License-Identifier: Apache-2.0
/*
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

#include <UnifiedVol.hpp>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>

using namespace uv;

int main(int argc, char* argv[])
{
    try
    {
        const std::filesystem::path path =
            (argc > 1) ? std::filesystem::path{argv[1]}
                       : std::filesystem::path{"data/VolSurface_SPY_04072011.csv"};

        initialize();

        using Real = double;

        core::MarketData<Real> marketData{
            .interestRate = 0.0,
            .dividendYield = 0.0,
            .spot = 485.77548
        };

        const core::MarketState<Real> marketState{io::load::marketState(path, marketData)
        };

        // io::report::volatility(marketState);

        // --------------  SVI calibration -------------

        const core::VolSurface<Real> sviVolSurface{models::svi::buildSurface(marketState)
        };

        // io::report::volatility(sviVolSurface);

        // --------------  Heston calibration --------------

        // TODO

        // Paramterize multithreading
        // Testing convergence on balck scholes would be cool
        // 2. Clean up pricer class
        // 4. Optimize pricer
        // Tanhsinh static_assert(N % 2 == 0, "TanHSinH<N>: N must be even for unroll-by-2
        // integration.");

        // Be smart, benchmark char funcito first with prof to see what is expensive!

        const core::VolSurface<Real> hestonVolSurface{
            models::heston::buildSurface<Real>(sviVolSurface, marketState.interestCurve)
        };

        utils::ScopedTimer timer{};

        // 20407.349381 ms

        std::size_t N{40};
        for (std::size_t i{0}; i < N; ++i)
        {
            const core::VolSurface<Real> hestonVolSurface{
                models::heston::buildSurface<Real>(
                    sviVolSurface,
                    marketState.interestCurve
                )
            };
        }

        // io::report::volatility(hestonVolSurface);

        return 0;

        // sviVolSurface.printBSCall();

        //// --------------  Local Vol calibration --------------

        // constexpr std::size_t nT{100};
        // constexpr std::size_t nX{10};
        // constexpr Real xBound{2.5};

        //// TODO paramterize
        // auto payoff = [](Real S, Real K) -> Real
        //{
        //     return (S > K) ? (S - K) : 0.0;
        // };

        // auto lvPricer =
        //     std::make_unique<localvol::Pricer<Real, nT, nX, PchipInterpolator<Real>>>(
        //         payoff,
        //         xBound
        //     );

        // Optimizer<Policy<void, ceres::LEVENBERG_MARQUARDT, ceres::DENSE_QR>>
        // lvOptimizer{
        //     opt::ceres::Config{

        //        .maxEval = 1000,
        //        .functionTol = 1e-14,
        //        .paramTol = 1e-14,
        //        .gradientTol = 1e-14,
        //        .verbose = false
        //    }
        //};

        // Surface<Real> localVolSurface{
        //     localvol::calibrator::calibrate<Real, nT, nX, PchipInterpolator<Real>>(
        //         sviVolSurface,
        //         *lvPricer,
        //         lvOptimizer,
        //         opt::WeightATM{.wATM = 8.0, .k0 = 0.3}
        //     )
        // };

        //// --------------  Heston calibration --------------

        // constexpr std::size_t HestonNodes{300};
        // const TanHSinH<Real, HestonNodes> quad{};

        // heston::Pricer<Real, HestonNodes> hestonPricer{
        //     std::make_shared<const TanHSinH<Real, HestonNodes>>(quad),
        //     {-2.0, 2.0}
        // };

        // Optimizer<Policy<void, ceres::LEVENBERG_MARQUARDT,
        // ceres::DENSE_NORMAL_CHOLESKY>>
        //     hestonOptimizer{opt::ceres::Config{

        //        .maxEval = 1000,
        //        .functionTol = 1e-12,
        //        .paramTol = 1e-12,
        //        .gradientTol = 1e-12,
        //        .paramNames = {"kappa", "theta", "sigma", "rho", "v0"},
        //        .verbose = false
        //    }};

        // hestonPricer.setParams(heston::calibrator::calibrate(
        //     sviVolSurface,
        //     hestonPricer,
        //     hestonOptimizer,
        //     opt::WeightATM{.wATM = 8.0, .k0 = 0.3}
        //));

        // const VolSurface<Real> hestonVolSurface{
        //     heston::calibrator::buildSurface<Real>(sviVolSurface, hestonPricer)
        // };

        // hestonVolSurface.printBSCall();

        // timer.LogTime<std::milli>();

        return EXIT_SUCCESS;
    }
    catch (const errors::UnifiedVolError& e)
    {
        std::cerr << e.what() << '\n';
        return 2;
    }
    catch (const std::exception& e)
    {
        std::cerr << "std::exception: " << e.what() << '\n';
        return 1;
    }
}