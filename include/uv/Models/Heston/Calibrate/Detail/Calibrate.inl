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

#include <ceres/ceres.h>

#include <algorithm>
#include <memory>
#include <span>

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM
)
{
    std::span<const T> maturities{volSurface.maturities()};

    return calibrate<T, N, Policy>(
        volSurface.maturities(),
        curve.interpolateDF(maturities),
        volSurface.forwards(),
        volSurface.strikes(),
        math::black::priceB76(volSurface, curve),
        pricer,
        optimizer,
        weightATM
    );
}

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& callM,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM
)
{
    detail::validateInputs<T>(maturities, discountFactors, forwards, strikes, callM);

    const Vector<double> maturitiesD{convertVector<double>(maturities)};
    const Vector<double> discountFactorsD{convertVector<double>(discountFactors)};
    const Vector<double> forwardsD{convertVector<double>(forwards)};
    const Vector<double> strikesD{convertVector<double>(strikes)};

    const auto callMD{callM.template as<double>()};

    const std::size_t numStrikes{strikesD.size()};

    optimizer
        .initialize(detail::initGuess(), detail::lowerBounds(), detail::upperBounds());

    optimizer.beginRun();

    Vector<double> bufferWeights(numStrikes);
    Vector<double> logKF(numStrikes);

    for (std::size_t i = 0; i < maturitiesD.size(); ++i)
    {

        std::span<const double> callMDRow{callMD[i]};
        const double t{maturitiesD[i]};
        const double dF{discountFactorsD[i]};
        const double F{forwardsD[i]};

        math::vol::logKF<double>(logKF, F, strikesD, true);

        opt::cost::weightsATM<double>(logKF, weightATM, bufferWeights);

        for (std::size_t j = 0; j < numStrikes; ++j)
        {
            optimizer.addResidualBlock(std::make_unique<detail::PriceResidualJac<T, N>>(
                t,
                dF,
                F,
                strikesD[j],
                callMDRow[j],
                bufferWeights[j],
                pricer
            ));
        }
    }

    std::span<const double> params{optimizer.solve()};

    return Params<T>{
        T(params[0]),
        T(params[1]),
        T(params[2]),
        T(params[3]),
        T(params[4])
    };
}

} // namespace uv::models::heston

namespace uv::models::heston::detail
{
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

template <std::floating_point T, std::size_t N>
struct PriceResidualJac final : public ceres::SizedCostFunction<1, 5>
{
    const double t_, dF_, F_, K_, callPriceMkt_, w_;
    const Pricer<T, N>* pricer_;

    PriceResidualJac(
        double t,
        double dF,
        double F,
        double K,
        double callPriceMkt,
        double w,
        const Pricer<T, N>& pricer
    ) noexcept
        : t_(t),
          dF_(dF),
          F_(F),
          K_(K),
          w_(w),
          callPriceMkt_(callPriceMkt),
          pricer_(&pricer)
    {
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)
        const override
    {
        const double* p = parameters[0];
        const auto pg =
            pricer_->callPriceWithGradient(p[0], p[1], p[2], p[3], p[4], t_, dF_, F_, K_);

        residuals[0] = double((pg[0] - callPriceMkt_) * w_);

        if ((jacobians != nullptr) && (jacobians[0] != nullptr))
        {
            double* J = jacobians[0];
            J[0] = double(pg[1] * w_);
            J[1] = double(pg[2] * w_);
            J[2] = double(pg[3] * w_);
            J[3] = double(pg[4] * w_);
            J[4] = double(pg[5] * w_);
        }
        return true;
    }
};
} // namespace uv::models::heston::detail