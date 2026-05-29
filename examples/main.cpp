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

#include "IO/Report.hpp"
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
        utils::ScopedTimer timer{};

        using Real = double;

        core::MarketData<Real> marketData{
            .interestRate = 0.0,
            .dividendYield = 0.0,
            .spot = 485.77548
        };

        // -------------- Market data -------------

        const core::MarketState<Real> marketState{io::load::marketState(path, marketData)
        };

        io::report::volatility(marketState);

        // --------------  SVI calibration -------------

        const core::VolSurface<Real> sviVolSurface{models::svi::buildSurface(marketState)
        };

        io::report::volatility(sviVolSurface);

        // --------------  Heston calibration --------------

        const core::VolSurface<Real> hestonVolSurface{
            models::heston::buildSurface<Real>(sviVolSurface, marketState.interestCurve)
        };

        io::report::volatility(hestonVolSurface);

        return 0;

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