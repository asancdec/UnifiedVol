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
#include "Base/Types.hpp"
#include "Math/Interpolation/Interpolator.hpp"

namespace uv::math::vol
{

template <std::floating_point T> T logKF(T F, T K, bool doValidate)
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(K);
        UV_REQUIRE_FINITE(F);

        UV_REQUIRE_NON_NEGATIVE(K);
        UV_REQUIRE_NON_NEGATIVE(F);
    }

    return std::log(K / F);
}

template <std::floating_point T>
void logKF(std::span<T> out, T F, std::span<const T> K, bool doValidate)
{
    if (doValidate)
    {
        UV_REQUIRE_SAME_SIZE(K, out);

        UV_REQUIRE_FINITE(K);
        UV_REQUIRE_FINITE(F);

        UV_REQUIRE_NON_NEGATIVE(K);
        UV_REQUIRE_NON_NEGATIVE(F);
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
        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(vol);

        UV_REQUIRE_NON_NEGATIVE(t);
        UV_REQUIRE_NON_NEGATIVE(vol);
    }

    return vol * vol * t;
}

template <std::floating_point T>
void totalVariance(std::span<T> out, T t, std::span<const T> vol, bool doValidate)
{
    if (doValidate)
    {
        UV_REQUIRE_SAME_SIZE(vol, out);

        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(vol);

        UV_REQUIRE_NON_NEGATIVE(t);
        UV_REQUIRE_NON_NEGATIVE(vol);
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

template <std::floating_point T>
void volFromTotalVariance(
    std::span<T> out,
    T t,
    std::span<const T> totalVariance,
    bool doValidate
)
{
    if (doValidate)
    {
        UV_REQUIRE_SAME_SIZE(totalVariance, out);

        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(totalVariance);

        UV_REQUIRE_NON_NEGATIVE(t);
        UV_REQUIRE_NON_NEGATIVE(totalVariance);
    }

    T invT{1.0 / t};

    for (std::size_t i{0}; i < totalVariance.size(); ++i)
    {
        out[i] = std::sqrt(totalVariance[i] * invT);
    }
}

template <std::floating_point T>
core::Matrix<T> volFromTotalVariance(
    const std::span<const T> t,
    const core::Matrix<T>& totalVariance,
    bool doValidate
)
{
    UV_REQUIRE_SAME_SIZE(t, totalVariance.rows());

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
        UV_REQUIRE_FINITE(vol);
        UV_REQUIRE_NON_NEGATIVE(vol);
    }

    return vol * vol;
}
template <std::floating_point T>
void variance(std::span<T> out, std::span<const T> vol, bool doValidate)
{
    if (doValidate)
    {
        UV_REQUIRE_SAME_SIZE(vol, out);
        UV_REQUIRE_FINITE(vol);
        UV_REQUIRE_NON_NEGATIVE(vol);
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
        UV_REQUIRE_NON_EMPTY(parameters);
        UV_REQUIRE_NON_EMPTY(logKF);

        UV_REQUIRE_FINITE(parameters);
        UV_REQUIRE_FINITE(logKF);

        UV_REQUIRE_STRICTLY_INCREASING(logKF);

        UV_REQUIRE_SAME_SIZE(logKF, parameters);
    }

    return math::interp::PchipInterpolator<T>{}(0.0, logKF, parameters);
}

template <std::floating_point T>
T impliedVol(T callPrice, T t, T dF, T F, T K, bool doValidate)
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(callPrice);
        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(dF);
        UV_REQUIRE_FINITE(F);
        UV_REQUIRE_FINITE(K);

        UV_REQUIRE_NON_NEGATIVE(callPrice);
        UV_REQUIRE_NON_NEGATIVE(t);
        UV_REQUIRE_NON_NEGATIVE(dF);
        UV_REQUIRE_NON_NEGATIVE(F);
        UV_REQUIRE_NON_NEGATIVE(K);
    }

    return static_cast<T>(detail::impliedVolJackelCall(
        static_cast<double>(callPrice),
        static_cast<double>(t),
        static_cast<double>(dF),
        static_cast<double>(F),
        static_cast<double>(K)
    ));
}

template <std::floating_point T>
void impliedVol(
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
        UV_REQUIRE_NON_EMPTY(callPrices);
        UV_REQUIRE_NON_EMPTY(strikes);

        UV_REQUIRE_SAME_SIZE(callPrices, strikes);
        UV_REQUIRE_SAME_SIZE(callPrices, out);
    }

    for (std::size_t i{0}; i < callPrices.size(); ++i)
    {
        out[i] = impliedVol(callPrices[i], t, dF, F, strikes[i], doValidate);
    }
}

template <std::floating_point T>
core::Matrix<T> impliedVol(
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
        UV_REQUIRE_NON_EMPTY(maturities);
        UV_REQUIRE_NON_EMPTY(discountFactors);
        UV_REQUIRE_NON_EMPTY(forwards);

        UV_REQUIRE_SAME_SIZE(maturities, discountFactors);
        UV_REQUIRE_SAME_SIZE(maturities, callPrices.rows());
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

template <std::floating_point T>
core::Matrix<T> impliedVol(
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
        volSurface.strikes()
    );
}

} // namespace uv::math::vol
