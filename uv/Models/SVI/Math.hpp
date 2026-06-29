// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>

namespace uv::models::svi
{

template <std::floating_point T> struct Params;

template <std::floating_point T>
T totalVariance(T a, T b, T rho, T m, T sigma, T k) noexcept;

template <std::floating_point T> T totalVariance(const Params<T>& params, T k) noexcept;

template <std::floating_point T> T gk(T a, T b, T rho, T m, T sigma, T k) noexcept;

template <std::floating_point T> T gk(const Params<T>& params, T k) noexcept;

} // namespace uv::models::svi

#include "Models/SVI/Detail/Math.inl"
