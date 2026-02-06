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
#include "Math/Interpolation/Interpolator.hpp"
#include "Math/Interpolation/Policies.hpp"
#include "Models/LocalVol/VarianceView.hpp"
#include "Utils/Aux/Errors.hpp"

#include <format>

#include <iostream>

namespace uv::models::localvol::calibrator
{

using core::Matrix;
using core::VolSurface;
using ErrorCode::InvalidArgument;

template <
    std::floating_point T,
    std::size_t NT,
    std::size_t NX,
    class Interpolator,
    typename Policy>
Surface<T> calibrate(
    const VolSurface<T>& volSurface,
    Pricer<T, NT, NX, Interpolator>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
)
{
    return calibrate<T, NT, NX, Interpolator>(
        volSurface.normalizedCallPrices(),
        volSurface.maturities(),
        volSurface.logKFMatrix(),
        volSurface.totVarMatrix(),
        pricer,
        optimizer,
        weightATM
    );
}

template <
    std::floating_point T,
    std::size_t NT,
    std::size_t NX,
    class Interpolator,
    typename Policy>
Surface<T> calibrate(
    const Matrix<T>& callM,
    const Vector<T>& maturities,
    const Matrix<T>& logKF,
    const Matrix<T>& totVar,
    Pricer<T, NT, NX, Interpolator>& pricer,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::WeightATM<double>& weightATM
)
{
    detail::validate(callM, maturities, logKF, totVar);

    std::size_t numMaturities{logKF.rows()};
    std::size_t numStrikes{logKF.cols()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        // Guess
        // detail::coldGuess<T>(maturities.front(), totVar[0], localVar);
    }

    // Vector<T> dydx(numStrikes);

    // typename Pricer<T, NT, NX, Interpolator>::derivatives_type derivPolicy{};
    // std::span<const T> logKFSlice{logKF[0]};
    // T tenor{maturities[0]};
    // derivPolicy(logKFSlice, localVar, dydx);
    // VarianceView<T> localVarView{logKFSlice, localVar, dydx};
    // Vector<T> prices{pricer.priceNormalized(tenor, logKFSlice, localVarView)};

    // for (const auto& p : callM[0])
    //{
    //     std::cout << p << ' ';
    // }
    // std::cout << "\n";

    return Surface<T>{maturities, logKF, logKF};
}

} // namespace uv::models::localvol::calibrator

namespace uv::models::localvol::calibrator::detail
{
template <std::floating_point T>
void validate(
    const Matrix<T>& callM,
    const Vector<T>& maturities,
    const Matrix<T>& logKF,
    const Matrix<T>& totVar
)
{

    const std::size_t numRows{callM.rows()};
    const std::size_t numCols{callM.cols()};

    UV_REQUIRE(
        numCols > 0 && numCols > 0,
        InvalidArgument,
        std::format("validate: callM must be non-empty, got {}x{}", numRows, numCols)
    );

    UV_REQUIRE(
        numRows == logKF.rows() && numCols == logKF.cols(),
        InvalidArgument,
        std::format(
            "validate: callM/logKF size mismatch — callM is "
            "{}x{}, logKF is {}x{}",
            numRows,
            numCols,
            logKF.rows(),
            logKF.cols()
        )
    );

    UV_REQUIRE(
        numRows == totVar.rows() && numCols == totVar.cols(),
        InvalidArgument,
        std::format(
            "validate: callM/totVar size mismatch — callM is "
            "{}x{}, totVar is {}x{}",
            numRows,
            numCols,
            totVar.rows(),
            totVar.cols()
        )
    );

    UV_REQUIRE(
        maturities.size() == numRows,
        InvalidArgument,
        std::format(
            "validate: row mismatch — maturities={}, rows={}",
            maturities.size(),
            numRows
        )
    );

    core::validateFinite<T>(maturities);

    for (std::size_t i = 0; i < numRows; ++i)
    {

        UV_REQUIRE(
            maturities[i] > 0.0,
            InvalidArgument,
            std::format("validate: maturities[{}]={} must be > 0", i, maturities[i])
        );

        if (i > 0)
        {
            UV_REQUIRE(
                maturities[i] > maturities[i - 1],
                InvalidArgument,
                std::format(
                    "validate: maturities must be strictly increasing, but "
                    "maturities[{}]={} <= maturities[{}]={}",
                    i,
                    maturities[i],
                    i - 1,
                    maturities[i - 1]
                )
            );
        }
    }

    for (std::size_t t = 0; t < numRows; ++t)
    {
        const std::span<const T> callMRow{callM[t]};
        const std::span<const T> logKFRow{logKF[t]};
        const std::span<const T> totVarRow{totVar[t]};

        core::validateFinite<T>(callMRow);
        core::validateFinite<T>(logKFRow);
        core::validateFinite<T>(totVarRow);

        for (std::size_t k = 0; k < numCols; ++k)
        {
            const T c{callMRow[k]};
            const T w{totVarRow[k]};

            UV_REQUIRE(
                c >= T(0),
                InvalidArgument,
                std::format("validate: callM({}, {})={} must be >= 0", t, k, c)
            );

            UV_REQUIRE(
                w >= T(0),
                InvalidArgument,
                std::format("validate: totVar({}, {})={} must be >= 0", t, k, w)
            );
        }

        for (std::size_t k = 1; k < numCols; ++k)
        {
            const T prev{logKFRow[k - 1]};
            const T curr{logKFRow[k]};

            UV_REQUIRE(
                curr > prev,
                InvalidArgument,
                std::format(
                    "validate: logKF row {} must be strictly increasing, but logKF({}, "
                    "{})={} <= logKF({}, {})={}",
                    t,
                    t,
                    k,
                    curr,
                    t,
                    k - 1,
                    prev
                )
            );
        }
    }

    for (std::size_t k = 0; k < numCols; ++k)
    {
        for (std::size_t t = 1; t < numRows; ++t)
        {
            const T prev{totVar[t - 1][k]};
            const T curr{totVar[t][k]};

            UV_REQUIRE(
                curr >= prev,
                InvalidArgument,
                std::format(
                    "validate: totVar must be non-decreasing in tenor at strike col {}, "
                    "but totVar({}, {})={} < totVar({}, {})={}",
                    k,
                    t,
                    k,
                    curr,
                    t - 1,
                    k,
                    prev
                )
            );
        }
    }
}

template <std::floating_point T>
Vector<double> coldGuess(T tenor, std::span<const T> totVar) noexcept
{
    std::size_t n{totVar.size()};
    T invTenor{1.0 / tenor};

    Vector<double> out(n);

    for (std::size_t i{0}; i < n; ++i)
    {
        out[i] = static_cast<double>(totVar[i] * invTenor);
    }

    return out;
}

} // namespace uv::models::localvol::calibrator::detail