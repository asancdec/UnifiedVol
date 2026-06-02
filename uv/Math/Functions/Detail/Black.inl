// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Base/Types.hpp"
#include "Math/Functions/Primitive.hpp"

#include <cmath>
#include <span>

namespace uv::math::black
{
template <std::floating_point T> core::Matrix<T>
priceB76(const core::VolSurface<T>& volSurface, const core::Curve<T>& curve, bool isCall)
{
    std::span<const T> t(volSurface.maturities());
    Vector<T> dF(curve.interpolateDF(t));
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

template <std::floating_point T> void priceB76(
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
    if (doValidate)
    {

        REQUIRE_SAME_SIZE(vol, K);
        REQUIRE_SAME_SIZE(vol, out);

        REQUIRE_FINITE(vol);
        REQUIRE_FINITE(t);
        REQUIRE_FINITE(dF);
        REQUIRE_FINITE(F);
        REQUIRE_FINITE(K);

        REQUIRE_POSITIVE(vol);
        REQUIRE_POSITIVE(t);
        REQUIRE_POSITIVE(dF);
        REQUIRE_POSITIVE(F);
        REQUIRE_POSITIVE(K);
    }

    for (std::size_t i{0}; i < vol.size(); ++i)
    {
        out[i] = priceB76(t, dF, F, vol[i], K[i], false, isCall);
    }
}

template <std::floating_point T>
T priceBS(T t, T r, T q, T vol, T S, T K, bool doValidate, bool isCall)
{
    if (doValidate)
    {
        REQUIRE_FINITE(vol);
        REQUIRE_FINITE(t);
        REQUIRE_FINITE(r);
        REQUIRE_FINITE(q);
        REQUIRE_FINITE(S);
        REQUIRE_FINITE(K);

        REQUIRE_POSITIVE(vol);
        REQUIRE_POSITIVE(S);
        REQUIRE_POSITIVE(K);
        REQUIRE_POSITIVE(t);
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
        REQUIRE_FINITE(vol);
        REQUIRE_FINITE(t);
        REQUIRE_FINITE(dF);
        REQUIRE_FINITE(F);
        REQUIRE_FINITE(K);

        REQUIRE_POSITIVE(vol);
        REQUIRE_POSITIVE(t);
        REQUIRE_POSITIVE(dF);
        REQUIRE_POSITIVE(F);
        REQUIRE_POSITIVE(K);
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

template <std::floating_point T> T vegaB76(T t, T dF, T F, T vol, T K) noexcept
{
    T d1{detail::d1FromForward(t, vol, F, K)};

    return dF * F * normalPDF(d1) * std::sqrt(t);
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
