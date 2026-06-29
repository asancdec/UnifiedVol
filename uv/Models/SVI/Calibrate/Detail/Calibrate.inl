// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Models/SVI/Calibrate/Detail/Constraints.hpp"
#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Models/SVI/Calibrate/Detail/Initialize.hpp"
#include "Models/SVI/Calibrate/Detail/Objective.hpp"

#include <cstddef>
#include <format>

namespace uv::models::svi::detail
{
template <std::floating_point T, opt::nlopt::Algorithm Algo> Params<T> calibrateSlice(
    T t,
    std::span<const double> logKF,
    std::span<const double> totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    const Params<T>* prevParams
);

template <std::floating_point T> void validateInputs(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance
);

template <std::floating_point T>
void logParams(const Params<T>& params, unsigned int valuePrec = 6);
} // namespace uv::models::svi::detail

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
            detail::logParams(sliceParams);

        surfaceParams.emplace_back(sliceParams);
    }

    return surfaceParams;
}

} // namespace uv::models::svi

namespace uv::models::svi::detail
{

template <std::floating_point T>
void logParams(const Params<T>& params, unsigned int valuePrec)
{
    INFO(std::format(
        "T={:.4f}, a={:.{}f}, b={:.{}f}, rho={:.{}f}, m={:.{}f}, sigma={:.{}f}",
        params.t,
        params.a,
        valuePrec,
        params.b,
        valuePrec,
        params.rho,
        valuePrec,
        params.m,
        valuePrec,
        params.sigma,
        valuePrec
    ));
}

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
    REQUIRE_NON_EMPTY(maturities);
    REQUIRE_FINITE(maturities);
    REQUIRE_NON_NEGATIVE(maturities);
    REQUIRE_STRICTLY_INCREASING(maturities);

    REQUIRE_SAME_SIZE(maturities, logKF.rows());
    REQUIRE_SAME_SIZE(maturities, totalVariance.rows());
    REQUIRE_SAME_SIZE(logKF.cols(), totalVariance.cols());

    for (std::size_t i{0}; i < maturities.size(); ++i)
    {

        std::span<const T> logKFSlice{logKF[i]};
        std::span<const T> totalVarianceSlice{totalVariance[i]};

        REQUIRE_NON_EMPTY(logKFSlice);
        REQUIRE_NON_EMPTY(totalVarianceSlice);

        REQUIRE_FINITE(logKFSlice);
        REQUIRE_FINITE(totalVarianceSlice);

        REQUIRE_STRICTLY_INCREASING(logKFSlice);

        REQUIRE_NON_NEGATIVE(totalVarianceSlice);
    }
}

} // namespace uv::models::svi::detail
