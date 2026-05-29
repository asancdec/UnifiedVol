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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#include "Base/Macros/Require.hpp"
#include "Core/Generate.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Math/LinearAlgebra/MatrixOps.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/NLoptAdapter.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <cmath>
#include <cstddef>
#include <span>

namespace uv::models::svi
{

template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::MarketState<T>& marketState, const Config& config)
{
    return buildSurface(marketState.volSurface, config);
}

template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Config& config)
{
    opt::nlopt::Optimizer<4, opt::nlopt::Algorithm::LD_SLSQP> nloptOptimizer{
        detail::makeNLoptConfig(config)
    };

    return buildSurface(
        volSurface,
        calibrate(volSurface, nloptOptimizer, config.printParams)
    );
}

template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Vector<Params<T>>& params)
{

    std::span<const T> maturities{volSurface.maturities()};

    UV_REQUIRE_SAME_SIZE(maturities, params);

    const core::Matrix<T> logKF{math::vol::logKF(volSurface)};

    return core::generateVolSurface<T>(
        volSurface,
        math::vol::volFromTotalVariance<T>(
            maturities,
            math::linear_algebra::generateIndexed<T>(
                volSurface.numMaturities(),
                volSurface.numStrikes(),
                [&](std::size_t i, std::size_t j)
                {
                    const Params<T>& p{params[i]};

                    return totalVariance(p.a, p.b, p.rho, p.m, p.sigma, logKF[i][j]);
                }
            )
        )
    );
}
} // namespace uv::models::svi