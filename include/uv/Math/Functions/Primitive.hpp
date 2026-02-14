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

#pragma once

#include "Base/Types.hpp"

#include <concepts>

namespace uv::math
{
template <std::floating_point T>
[[gnu::hot]] Complex<T> log1pComplex(Complex<T> z) noexcept;

template <std::floating_point T> [[gnu::hot]] T cosm1(T b) noexcept;

template <std::floating_point T>
[[gnu::hot]] Complex<T> expm1Complex(Complex<T> z) noexcept;

template <std::floating_point T> T normalCDF(T x) noexcept;

template <std::floating_point T> T normalPDF(T x) noexcept;
} // namespace uv::math

#include "Math/Functions/Detail/Primitive.inl"