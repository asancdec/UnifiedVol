// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"

#include <concepts>

namespace uv::math::black
{

template <std::floating_point T> core::Matrix<T> priceB76(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    bool isCall = true
);

template <std::floating_point T> void priceB76(
    std::span<T> out,
    T t,
    T dF,
    T F,
    std::span<const T> vol,
    std::span<const T> K,
    bool doValidate = true,
    bool isCall = true
);

template <std::floating_point T>
T priceB76(T t, T dF, T F, T vol, T K, bool doValidate = true, bool isCall = true);

template <std::floating_point T>
T priceBS(T t, T r, T q, T vol, T S, T K, bool doValidate = true, bool isCall = true);

template <std::floating_point T> T vegaB76(T t, T dF, T F, T vol, T K) noexcept;

} // namespace uv::math::black

#include "Math/Functions/Detail/Black.inl"
