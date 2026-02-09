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

#include "Base/Macros/Require.hpp"
#include "Core/Matrix.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Math/LinearAlgebra/MatrixOps.hpp"
#include "Models/Heston/Calibrate/CeresAdapter.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Calibrate/Detail/ResidualCost.hpp"

#include <algorithm>
#include <span>

namespace uv::models::heston::calibrate
{

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config
)
{
    price::Pricer<T, N> pricer{};
    opt::ceres::Optimizer<Policy> optimizer{detail::makeOptimizer(config)};

    return detail::calibrate<T, N, Mode, Policy>(
        volSurface,
        curve,
        optimizer,
        config.weightATM,
        pricer
    );
}

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config,
    price::Pricer<T, N>& pricer
)
{
    opt::ceres::Optimizer<Policy> optimizer{detail::makeOptimizer(config)};

    return detail::calibrate<T, N, Mode, Policy>(
        volSurface,
        curve,
        optimizer,
        config.weightATM,
        pricer
    );
}

} // namespace uv::models::heston::calibrate

namespace uv::models::heston::calibrate::detail
{

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
)
{
    std::span<const T> maturities{volSurface.maturities()};

    return calibrate<T, N, Mode, Policy>(
        volSurface.maturities(),
        curve.interpolateDF(maturities),
        volSurface.forwards(),
        volSurface.strikes(),
        math::black::priceB76(volSurface, curve),
        optimizer,
        weightATM,
        pricer
    );
}

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
)
{
    if constexpr (std::is_same_v<T, double>)
    {
        return calibrateDouble<T, N, Mode, Policy>(
            maturities,
            discountFactors,
            forwards,
            strikes,
            callM,
            optimizer,
            weightATM,
            pricer
        );
    }

    return calibrateDouble<T, N, Mode, Policy>(
               convertVector<double>(maturities),
               convertVector<double>(discountFactors),
               convertVector<double>(forwards),
               convertVector<double>(strikes),
               callM.template as<double>(),
               optimizer,
               weightATM,
               pricer
    )
        .template as<T>();
}

template <
    std::floating_point CalcT,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<double> calibrateDouble(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& callM,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
)
{
    validateInputs(maturities, discountFactors, forwards, strikes, callM);

    const std::size_t numStrikes{strikes.size()};

    optimizer
        .initialize(detail::initGuess(), detail::lowerBounds(), detail::upperBounds());

    optimizer.beginRun();

    Vector<double> bufferWeights(numStrikes);
    Vector<double> logKF(numStrikes);

    for (std::size_t i = 0; i < maturities.size(); ++i)
    {

        std::span<const double> callMRow{callM[i]};
        const double t{maturities[i]};
        const double dF{discountFactors[i]};
        const double F{forwards[i]};

        math::vol::logKF<double>(logKF, F, strikes, true);

        opt::cost::weightsATM<double>(logKF, weightATM, bufferWeights);

        for (std::size_t j = 0; j < numStrikes; ++j)
        {
            optimizer.addResidualBlock(makeCost<Mode, CalcT, N>(
                t,
                dF,
                F,
                strikes[j],
                callMRow[j],
                bufferWeights[j],
                pricer
            ));
        }
    }

    return Params<double>{optimizer.solve()};
}

template <std::floating_point T>
void validateInputs(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM
)
{
    UV_REQUIRE_NON_EMPTY(maturities);
    UV_REQUIRE_NON_EMPTY(discountFactors);
    UV_REQUIRE_NON_EMPTY(forwards);
    UV_REQUIRE_NON_EMPTY(strikes);

    UV_REQUIRE_FINITE(maturities);
    UV_REQUIRE_FINITE(discountFactors);
    UV_REQUIRE_FINITE(forwards);
    UV_REQUIRE_FINITE(strikes);

    UV_REQUIRE_POSITIVE(maturities);
    UV_REQUIRE_POSITIVE(discountFactors);

    UV_REQUIRE_SAME_SIZE(maturities, forwards);
    UV_REQUIRE_SAME_SIZE(maturities, discountFactors);
    UV_REQUIRE_SAME_SIZE(maturities, callM.rows());
    UV_REQUIRE_SAME_SIZE(strikes, callM.cols());

    for (std::size_t i{0}; i < maturities.size(); ++i)
    {
        std::span<const T> callMRow{callM[i]};

        UV_REQUIRE_NON_EMPTY(callMRow);
        UV_REQUIRE_FINITE(callMRow);
        UV_REQUIRE_POSITIVE(callMRow);
    }
}
} // namespace uv::models::heston::calibrate::detail