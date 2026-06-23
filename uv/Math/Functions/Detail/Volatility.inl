// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Base/Types.hpp"
#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Interpolation/Hermite/Interpolator.hpp"
#include "Math/LinearAlgebra/VectorOps.hpp"

#include <cmath>
#include <cstddef>
#include <span>

namespace uv::math::vol::detail
{
double impliedVolJackelCall(double callPrice, double t, double dF, double F, double K);
} // namespace uv::math::vol::detail

namespace uv::math::vol
{

template <std::floating_point T> T logKF(T F, T K, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_FINITE(K);
        REQUIRE_FINITE(F);

        REQUIRE_POSITIVE(K);
        REQUIRE_POSITIVE(F);
    }

    return std::log(K / F);
}

template <std::floating_point T>
void logKF(std::span<T> out, T F, std::span<const T> K, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_SAME_SIZE(K, out);

        REQUIRE_FINITE(K);
        REQUIRE_FINITE(F);

        REQUIRE_POSITIVE(K);
        REQUIRE_POSITIVE(F);
    }

    for (std::size_t i{0}; i < K.size(); ++i)
    {
        out[i] = logKF(F, K[i], false);
    }
}

template <std::floating_point T>
core::Matrix<T> logKF(const core::VolSurface<T>& volSurface, bool doValidate)
{
    std::size_t numMaturities{volSurface.numMaturities()};

    std::span<const T> forwards(volSurface.forwards());
    std::span<const T> strikes{volSurface.strikes()};

    core::Matrix<T> out{numMaturities, volSurface.numStrikes()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        logKF(out[i], forwards[i], strikes, doValidate);
    }

    return out;
}

template <std::floating_point T> T totalVariance(T t, T vol, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_FINITE(t);
        REQUIRE_FINITE(vol);

        REQUIRE_NON_NEGATIVE(t);
        REQUIRE_NON_NEGATIVE(vol);
    }

    return vol * vol * t;
}

template <std::floating_point T>
void totalVariance(std::span<T> out, T t, std::span<const T> vol, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_SAME_SIZE(vol, out);

        REQUIRE_FINITE(t);
        REQUIRE_FINITE(vol);

        REQUIRE_NON_NEGATIVE(t);
        REQUIRE_NON_NEGATIVE(vol);
    }

    for (std::size_t i{0}; i < vol.size(); ++i)
    {
        out[i] = totalVariance(t, vol[i], false);
    }
}

template <std::floating_point T>
core::Matrix<T> totalVariance(const core::VolSurface<T>& volSurface, bool doValidate)
{
    std::size_t numMaturities{volSurface.numMaturities()};

    std::span<const T> maturities(volSurface.maturities());
    const core::Matrix<T>& vol{volSurface.vol()};

    core::Matrix<T> out{numMaturities, volSurface.numStrikes()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        totalVariance(out[i], maturities[i], vol[i], doValidate);
    }

    return out;
}

template <std::floating_point T> void volFromTotalVariance(
    std::span<T> out,
    T t,
    std::span<const T> totalVariance,
    bool doValidate
)
{
    if (doValidate)
    {
        REQUIRE_SAME_SIZE(totalVariance, out);

        REQUIRE_FINITE(t);
        REQUIRE_FINITE(totalVariance);

        REQUIRE_POSITIVE(t);
        REQUIRE_NON_NEGATIVE(totalVariance);
    }

    const T invT{1.0 / t};
    const Vector<T> variancePerUnitTime{
        math::linear_algebra::multiply<T>(totalVariance, invT)
    };

    math::linear_algebra::squareRootInplace<T>(out, variancePerUnitTime);
}

template <std::floating_point T> core::Matrix<T> volFromTotalVariance(
    const std::span<const T> t,
    const core::Matrix<T>& totalVariance,
    bool doValidate
)
{
    REQUIRE_SAME_SIZE(t, totalVariance.rows());

    const std::size_t n{t.size()};

    core::Matrix<T> out(n, totalVariance.cols());

    for (std::size_t i{0}; i < n; ++i)
    {
        volFromTotalVariance(out[i], t[i], totalVariance[i], doValidate);
    }

    return out;
}

template <std::floating_point T> T variance(T vol, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_FINITE(vol);
        REQUIRE_NON_NEGATIVE(vol);
    }

    return vol * vol;
}
template <std::floating_point T>
void variance(std::span<T> out, std::span<const T> vol, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_SAME_SIZE(vol, out);
        REQUIRE_FINITE(vol);
        REQUIRE_NON_NEGATIVE(vol);
    }

    for (std::size_t i{0}; i < vol.size(); ++i)
    {
        out[i] = variance(vol[i], false);
    }
}

template <std::floating_point T>
core::Matrix<T> variance(const core::VolSurface<T>& volSurface, bool doValidate)
{
    std::size_t numMaturities{volSurface.numMaturities()};
    const core::Matrix<T>& vol{volSurface.vol()};

    core::Matrix<T> out{numMaturities, volSurface.numStrikes()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        variance(out[i], vol[i], doValidate);
    }

    return out;
}

template <std::floating_point T>
T atmParameter(std::span<const T> parameters, std::span<const T> logKF, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_NON_EMPTY(parameters);
        REQUIRE_NON_EMPTY(logKF);

        REQUIRE_FINITE(parameters);
        REQUIRE_FINITE(logKF);

        REQUIRE_STRICTLY_INCREASING(logKF);

        REQUIRE_SAME_SIZE(logKF, parameters);
    }

    return math::interp::hermite::PchipInterpolator<T>{}(0.0, logKF, parameters);
}

template <std::floating_point T>
T impliedVol(T callPrice, T t, T dF, T F, T K, bool doValidate)
{
    if (doValidate)
    {
        REQUIRE_FINITE(callPrice);
        REQUIRE_FINITE(t);
        REQUIRE_FINITE(dF);
        REQUIRE_FINITE(F);
        REQUIRE_FINITE(K);

        REQUIRE_NON_NEGATIVE(callPrice);
        REQUIRE_NON_NEGATIVE(t);
        REQUIRE_NON_NEGATIVE(dF);
        REQUIRE_NON_NEGATIVE(F);
        REQUIRE_NON_NEGATIVE(K);
    }

    return static_cast<T>(detail::impliedVolJackelCall(
        static_cast<double>(callPrice),
        static_cast<double>(t),
        static_cast<double>(dF),
        static_cast<double>(F),
        static_cast<double>(K)
    ));
}

template <std::floating_point T> void impliedVol(
    std::span<T> out,
    std::span<const T> callPrices,
    T t,
    T dF,
    T F,
    std::span<const T> strikes,
    bool doValidate
)
{

    if (doValidate)
    {
        REQUIRE_NON_EMPTY(callPrices);
        REQUIRE_NON_EMPTY(strikes);

        REQUIRE_SAME_SIZE(callPrices, strikes);
        REQUIRE_SAME_SIZE(callPrices, out);
    }

    for (std::size_t i{0}; i < callPrices.size(); ++i)
    {
        out[i] = impliedVol(callPrices[i], t, dF, F, strikes[i], doValidate);
    }
}

template <std::floating_point T> core::Matrix<T> impliedVol(
    const core::Matrix<T>& callPrices,
    std::span<const T> maturities,
    std::span<const T> discountFactors,
    std::span<const T> forwards,
    std::span<const T> strikes,
    bool doValidate
)
{

    if (doValidate)
    {
        REQUIRE_NON_EMPTY(maturities);
        REQUIRE_NON_EMPTY(discountFactors);
        REQUIRE_NON_EMPTY(forwards);

        REQUIRE_SAME_SIZE(maturities, discountFactors);
        REQUIRE_SAME_SIZE(maturities, callPrices.rows());
    }

    std::size_t numMaturities{maturities.size()};

    core::Matrix<T> out{numMaturities, strikes.size()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        impliedVol(
            out[i],
            callPrices[i],
            maturities[i],
            discountFactors[i],
            forwards[i],
            strikes,
            doValidate
        );
    }

    return out;
}

template <std::floating_point T> core::Matrix<T> impliedVol(
    const core::Matrix<T>& callPrices,
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    bool doValidate
)
{

    std::span<const T> maturities{volSurface.maturities()};
    const Vector<T> discountFactors{curve.interpolateDF(maturities)};

    return impliedVol<T>(
        callPrices,
        maturities,
        discountFactors,
        volSurface.forwards(),
        volSurface.strikes(),
        doValidate
    );
}

} // namespace uv::math::vol
