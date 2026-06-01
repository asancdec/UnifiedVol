// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <span>

namespace uv::models::heston
{

template <std::floating_point T> struct Params
{
    T kappa;
    T theta;
    T sigma;
    T rho;
    T v0;

    constexpr Params(T kappa_, T theta_, T sigma_, T rho_, T v0_) noexcept;

    explicit Params(std::span<const double> params);

    template <std::floating_point U> constexpr Params<U> as() const noexcept;
};

} // namespace uv::models::heston

#include "Models/Heston/Detail/Params.inl"