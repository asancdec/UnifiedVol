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
#include "IO/Report.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Models/SVI/Calibrate/Detail/Constraints.hpp"
#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Models/SVI/Calibrate/Detail/Initialize.hpp"
#include "Models/SVI/Calibrate/Detail/Objective.hpp"

#include <cstddef>

namespace uv::models::svi
{
template <std::floating_point T, opt::nlopt::Algorithm Algo> Vector<Params<T>> calibrate(
    const core::VolSurface<T>& volSurface,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    bool printParams
)
{
    return calibrate<T, Algo>(
        volSurface.maturities(),
        math::vol::logKF(volSurface),
        math::vol::totalVariance(volSurface),
        prototype,
        printParams
    );
}

template <std::floating_point T, opt::nlopt::Algorithm Algo> Vector<Params<T>> calibrate(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    bool printParams
)
{
    detail::validateInputs<T>(maturities, logKF, totalVariance);

    const auto logKFD{logKF.template as<double>()};
    const auto totalVarianceD{totalVariance.template as<double>()};

    const std::size_t numMaturities{maturities.size()};

    Vector<Params<T>> surfaceParams;
    surfaceParams.reserve(numMaturities);

    for (std::size_t i = 0; i < numMaturities; ++i)
    {

        Params<T> sliceParams{detail::calibrateSlice<T>(
            maturities[i],
            logKFD[i],
            totalVarianceD[i],
            prototype,
            (i == 0) ? nullptr : &surfaceParams.back()
        )};

        if (printParams)
            io::report::sviParams(sliceParams);

        surfaceParams.emplace_back(sliceParams);
    }

    return surfaceParams;
}

} // namespace uv::models::svi

namespace uv::models::svi::detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo> Params<T> calibrateSlice(
    const T t,
    std::span<const double> logKF,
    std::span<const double> totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    const Params<T>* prevParams
)
{
    const SliceData sliceData(logKF, totalVariance);

    const double atmTotalVariance{sliceData.atmTotalVariance};

    opt::nlopt::Optimizer optimizer{prototype.fresh()};
    optimizer.setUserValue(atmTotalVariance);

    setGuessBounds(optimizer, prevParams, sliceData);

    SliceConstraints c;
    addCalendarConstraints(optimizer, c, prevParams, logKF, sliceData);

    addWMinConstraint(optimizer);
    addMinSlopeConstraint(optimizer);
    addMaxSlopeConstraint(optimizer);

    ConvexityMContext convexityCtx;
    addConvexityConstraints(optimizer, convexityCtx, logKF, atmTotalVariance);

    ObjectiveContexts obj{logKF, totalVariance, atmTotalVariance};

    setMinObjective(optimizer, obj);

    return Params<T>{t, optimizer.optimize(), atmTotalVariance};
}

template <std::floating_point T> void validateInputs(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance
)
{
    UV_REQUIRE_NON_EMPTY(maturities);
    UV_REQUIRE_FINITE(maturities);
    UV_REQUIRE_NON_NEGATIVE(maturities);
    UV_REQUIRE_STRICTLY_INCREASING(maturities);

    UV_REQUIRE_SAME_SIZE(maturities, logKF.rows());
    UV_REQUIRE_SAME_SIZE(maturities, totalVariance.rows());
    UV_REQUIRE_SAME_SIZE(logKF.cols(), totalVariance.cols());

    for (std::size_t i{0}; i < maturities.size(); ++i)
    {

        std::span<const T> logKFSlice{logKF[i]};
        std::span<const T> totalVarianceSlice{totalVariance[i]};

        UV_REQUIRE_NON_EMPTY(logKFSlice);
        UV_REQUIRE_NON_EMPTY(totalVarianceSlice);

        UV_REQUIRE_FINITE(logKFSlice);
        UV_REQUIRE_FINITE(totalVarianceSlice);

        UV_REQUIRE_NON_NEGATIVE(totalVarianceSlice);
    }
}

} // namespace uv::models::svi::detail
