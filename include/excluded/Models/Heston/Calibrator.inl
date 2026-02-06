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

#include "Core/Functions.hpp"
#include "Core/Matrix/Functions.hpp"

#include <ceres/ceres.h>

#include <algorithm>
#include <memory>
#include <span>

namespace uv::models::heston::calibrator
{

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
)
{
    return calibrate<T, N, Policy>(
        volSurface.maturities(),
        volSurface.strikes(),
        volSurface.forwards(),
        volSurface.rates(),
        volSurface.callPrices(),
        pricer,
        optimizer,
        weightATM
    );
}

template <std::floating_point T, std::size_t N, typename Policy>
Params<T> calibrate(
    const Vector<T>& maturities,
    const Vector<T>& strikes,
    const Vector<T>& forwards,
    const Vector<T>& rates,
    const core::Matrix<T>& callM,
    Pricer<T, N>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
)
{

    detail::validateInputs<T>(maturities, strikes, forwards, rates, callM);

    const Vector<double> maturitiesD{core::convertVector<double>(maturities)};
    const Vector<double> strikesD{core::convertVector<double>(strikes)};
    const Vector<double> forwardsD{core::convertVector<double>(forwards)};
    const Vector<double> ratesD{core::convertVector<double>(rates)};
    const core::Matrix<double> callMD{core::convertMatrix<double>(callM)};

    const std::size_t numMaturities{maturitiesD.size()};
    const std::size_t numStrikes{strikesD.size()};

    optimizer
        .initialize(detail::initGuess(), detail::lowerBounds(), detail::upperBounds());

    optimizer.beginRun();

    Vector<double> bufferWeights(numStrikes);
    Vector<double> logKF(numStrikes);

    for (std::size_t i = 0; i < numMaturities; ++i)
    {

        std::span<const double> callMDRow{callMD[i]};
        const double t{maturitiesD[i]};
        const double F{forwardsD[i]};
        const double r{ratesD[i]};

        for (std::size_t j{0}; j < numStrikes; ++j)
            logKF[j] = std::log(strikesD[j] / F);

        opt::weightsATM<double>(logKF, weightATM, bufferWeights);

        for (std::size_t j = 0; j < numStrikes; ++j)
        {
            optimizer.addResidualBlock(std::make_unique<detail::PriceResidualJac<T, N>>(
                t,
                F,
                r,
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

template <std::floating_point T, std::size_t N>
core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Pricer<T, N>& pricer)
{

    const std::size_t numMaturities{volSurface.numMaturities()};
    const std::size_t numStrikes{volSurface.numStrikes()};
    const Vector<T>& maturities{volSurface.maturities()};
    const Vector<T>& strikes{volSurface.strikes()};
    const Vector<T>& forwards{volSurface.forwards()};
    const Vector<T>& rates{volSurface.rates()};

    core::VolSurface<T> hestonVolSurface{volSurface};
    hestonVolSurface.setCallPrices(core::generateIndexed<T>(
        numMaturities,
        numStrikes,
        [&](std::size_t i, std::size_t j)
        {
            return pricer.callPrice(maturities[i], forwards[i], rates[i], strikes[j]);
        }
    ));

    return hestonVolSurface;
}

} // namespace uv::models::heston::calibrator

namespace uv::models::heston::calibrator::detail
{
template <std::floating_point T>
void validateInputs(
    const Vector<T>& maturities,
    const Vector<T>& strikes,
    const Vector<T>& forwards,
    const Vector<T>& rates,
    const core::Matrix<T>& callM
)
{
    UV_REQUIRE(
        !maturities.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: maturities is empty"
    );

    UV_REQUIRE(
        !strikes.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: strikes is empty"
    );

    UV_REQUIRE(
        !forwards.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: forwards is empty"
    );

    UV_REQUIRE(
        !rates.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: rates is empty"
    );

    UV_REQUIRE(
        !callM.empty(),
        ErrorCode::InvalidArgument,
        "validateInputs: callM is empty"
    );

    const std::size_t numMaturities{maturities.size()};
    const std::size_t numStrikes{strikes.size()};

    UV_REQUIRE(
        forwards.size() == numMaturities,
        ErrorCode::InvalidArgument,
        "validateInputs: forwards size must equal maturities size"
    );

    UV_REQUIRE(
        rates.size() == numMaturities,
        ErrorCode::InvalidArgument,
        "validateInputs: rates size must equal maturities size"
    );

    UV_REQUIRE(
        callM.rows() == numMaturities,
        ErrorCode::InvalidArgument,
        "validateInputs: callM rows must equal number of maturities"
    );

    UV_REQUIRE(
        callM.cols() == numStrikes,
        ErrorCode::InvalidArgument,
        "validateInputs: callM columns must equal strikes size"
    );

    for (std::size_t i = 0; i < numMaturities; ++i)
    {
        UV_REQUIRE(
            maturities[i] > Real(0),
            ErrorCode::InvalidArgument,
            "validateInputs: maturities must be > 0"
        );

        if (i > 0)
        {
            UV_REQUIRE(
                maturities[i] > maturities[i - 1],
                ErrorCode::InvalidArgument,
                "validateInputs: maturities must be strictly increasing"
            );
        }
    }

    for (std::size_t j = 0; j < numStrikes; ++j)
    {
        UV_REQUIRE(
            strikes[j] > Real(0),
            ErrorCode::InvalidArgument,
            "validateInputs: strikes must be > 0"
        );

        if (j > 0)
        {
            UV_REQUIRE(
                strikes[j] > strikes[j - 1],
                ErrorCode::InvalidArgument,
                "validateInputs: strikes must be strictly increasing"
            );
        }
    }

    for (std::size_t i = 0; i < numMaturities; ++i)
    {
        UV_REQUIRE(
            forwards[i] > Real(0),
            ErrorCode::InvalidArgument,
            "validateInputs: forwards must be > 0"
        );
    }
}

template <std::floating_point T, std::size_t N>
struct PriceResidualJac final : public ceres::SizedCostFunction<1, 5>
{
    const double T_, F_, r_, K_, callPriceMkt_, w_;
    const Pricer<T, N>* pricer_;

    PriceResidualJac(
        double t,
        double F,
        double r,
        double K,
        double callPriceMkt,
        double w,
        const Pricer<T, N>& pricer
    ) noexcept
        : T_(t),
          F_(F),
          r_(r),
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
            pricer_->callPriceWithGradient(p[0], p[1], p[2], p[3], p[4], T_, F_, r_, K_);

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
} // namespace uv::models::heston::calibrator::detail