// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>

namespace uv::math
{

template <std::floating_point T> [[gnu::hot]] Complex<T> invComplex(Complex<T>) noexcept;

template <std::floating_point T>
[[gnu::hot]] Complex<T> log1pComplex(Complex<T>) noexcept;

template <std::floating_point T> [[gnu::hot]] T cosm1(T) noexcept;

template <std::floating_point T>
[[gnu::hot]] Complex<T> expm1Complex(Complex<T>) noexcept;

template <std::floating_point T> T normalCDF(T) noexcept;

template <std::floating_point T> T normalPDF(T) noexcept;
} // namespace uv::math

#include "Math/Functions/Detail/Primitive.inl"
