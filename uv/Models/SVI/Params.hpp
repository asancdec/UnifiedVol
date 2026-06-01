// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <span>

namespace uv::models::svi
{

template <std::floating_point T> struct Params
{
    T t;

    T a;
    T b;
    T rho;
    T m;
    T sigma;

    Params(T t_, T a_, T b_, T rho_, T m_, T sigma_) noexcept;

    Params(T t_, std::span<const double> params, double atmTotalVariance) noexcept;

    template <std::floating_point U> Params<U> as() const noexcept;
};
} // namespace uv::models::svi

#include "Models/SVI/Detail/Params.inl"