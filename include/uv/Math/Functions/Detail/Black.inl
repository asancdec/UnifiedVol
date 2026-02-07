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

#include "Base/Alias.hpp"
#include "Base/Macros/Require.hpp"
#include "Math/Functions/Primitive.hpp"

#include <cmath>
#include <numbers>
#include <span>

namespace uv::math::black
{
template <std::floating_point T>
core::Matrix<T>
priceB76(const core::VolSurface<T>& volSurface, const core::Curve<T>& curve, bool isCall)
{
    std::span<const T> t(volSurface.maturities());
    Vector<T> dF(curve.discountFactor(t));
    std::span<const T> F(volSurface.forwards());
    std::span<const T> K(volSurface.strikes());
    const core::Matrix<T>& vol{volSurface.vol()};

    std::size_t numMaturities{volSurface.numMaturities()};
    core::Matrix<T> out{numMaturities, volSurface.numStrikes()};

    for (std::size_t i{0}; i < numMaturities; ++i)
    {
        priceB76(out[i], t[i], dF[i], F[i], vol[i], K, true, isCall);
    }

    return out;
}

template <std::floating_point T>
void priceB76(
    std::span<T> out,
    T t,
    T dF,
    T F,
    std::span<const T> vol,
    std::span<const T> K,
    bool doValidate,
    bool isCall
)
{
    std::size_t n{vol.size()};

    if (doValidate)
    {

        UV_REQUIRE_SAME_SIZE(n, K.size());
        UV_REQUIRE_SAME_SIZE(n, out.size());

        UV_REQUIRE_FINITE(vol);
        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(dF);
        UV_REQUIRE_FINITE(F);
        UV_REQUIRE_FINITE(K);

        UV_REQUIRE_POSITIVE(vol);
        UV_REQUIRE_POSITIVE(t);
        UV_REQUIRE_POSITIVE(dF);
        UV_REQUIRE_POSITIVE(F);
        UV_REQUIRE_POSITIVE(K);
    }

    for (std::size_t i{0}; i < n; ++i)
    {
        out[i] = priceB76(t, dF, F, vol[i], K[i], false, isCall);
    }
}

template <std::floating_point T>
T priceBS(T t, T r, T q, T vol, T S, T K, bool doValidate, bool isCall)
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(vol);
        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(r);
        UV_REQUIRE_FINITE(q);
        UV_REQUIRE_FINITE(S);
        UV_REQUIRE_FINITE(K);

        UV_REQUIRE_POSITIVE(vol);
        UV_REQUIRE_POSITIVE(S);
        UV_REQUIRE_POSITIVE(K);
        UV_REQUIRE_POSITIVE(t);
    }

    const T d1{detail::d1(t, r, q, vol, S, K)};
    const T d2{detail::d2(vol, t, d1)};

    if (isCall)
    {
        return S * std::exp(-q * t) * normalCDF(d1) -
               K * std::exp(-r * t) * normalCDF(d2);
    }
    else
    {
        return K * std::exp(-r * t) * normalCDF(-d2) -
               S * std::exp(-q * t) * normalCDF(-d1);
    }
}

template <std::floating_point T>
T priceB76(T t, T dF, T F, T vol, T K, bool doValidate, bool isCall)
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(vol);
        UV_REQUIRE_FINITE(t);
        UV_REQUIRE_FINITE(dF);
        UV_REQUIRE_FINITE(F);
        UV_REQUIRE_FINITE(K);

        UV_REQUIRE_POSITIVE(vol);
        UV_REQUIRE_POSITIVE(t);
        UV_REQUIRE_POSITIVE(dF);
        UV_REQUIRE_POSITIVE(F);
        UV_REQUIRE_POSITIVE(K);
    }

    const T d1{detail::d1FromForward(t, vol, F, K)};
    const T d2{detail::d2(vol, t, d1)};

    if (isCall)
    {
        return dF * (F * normalCDF(d1) - K * normalCDF(d2));
    }
    else
    {
        return dF * (K * normalCDF(-d2) - F * normalCDF(-d1));
    }
}

} // namespace uv::math::black

namespace uv::math::black::detail
{

template <std::floating_point T> T d1(T t, T r, T q, T vol, T S, T K) noexcept
{
    return std::fma(t, (r - q + vol * vol * 0.5), std::log(S / K)) / (std::sqrt(t) * vol);
}

template <std::floating_point T> T d2(T vol, T t, T d1) noexcept
{
    return std::fma(-vol, std::sqrt(t), d1);
}

template <std::floating_point T> T d1FromForward(T t, T vol, T F, T K) noexcept
{
    return std::fma(t, (vol * vol * 0.5), std::log(F / K)) / (std::sqrt(t) * vol);
}

} // namespace uv::math::black::detail
