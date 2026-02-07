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

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"

#include <concepts>

namespace uv::math::black
{

template <std::floating_point T>
core::Matrix<T> priceB76(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    bool isCall = true
);

template <std::floating_point T>
void priceB76(
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

namespace detail
{

template <std::floating_point T> T d1(T t, T r, T q, T vol, T S, T K) noexcept;

template <std::floating_point T> T d1FromForward(T t, T vol, T F, T K) noexcept;

template <std::floating_point T> T d2(T vol, T t, T d1) noexcept;

} // namespace detail

} // namespace uv::math::black

#include "Math/Functions/Detail/Black.inl"