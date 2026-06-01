// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>

namespace uv::math
{

template <std::floating_point T>
[[gnu::hot]] Complex<T> invComplex(Complex<T> z) noexcept;

template <std::floating_point T>
[[gnu::hot]] Complex<T> log1pComplex(Complex<T> z) noexcept;

template <std::floating_point T> [[gnu::hot]] T cosm1(T b) noexcept;

template <std::floating_point T>
[[gnu::hot]] Complex<T> expm1Complex(Complex<T> z) noexcept;

template <std::floating_point T> T normalCDF(T x) noexcept;

template <std::floating_point T> T normalPDF(T x) noexcept;
} // namespace uv::math

#include "Math/Functions/Detail/Primitive.inl"