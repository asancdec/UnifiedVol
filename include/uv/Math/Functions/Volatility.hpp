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

#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"

#include <concepts>
#include <span>

namespace uv::math::vol
{

template <std::floating_point T> T logKF(T F, T K, bool doValidate = true);

template <std::floating_point T>
void logKF(std::span<T> out, T F, std::span<const T> K, bool doValidate = true);

template <std::floating_point T>
core::Matrix<T> logKF(const core::VolSurface<T>& volSurface, bool doValidate = true);

template <std::floating_point T> T totalVariance(T t, T vol, bool doValidate);

template <std::floating_point T>
void totalVariance(std::span<T> out, T t, std::span<const T> vol, bool doValidate = true);

template <std::floating_point T>
core::Matrix<T>
totalVariance(const core::VolSurface<T>& volSurface, bool doValidate = true);

template <std::floating_point T>
void volFromTotalVariance(
    std::span<T> out,
    T t,
    std::span<const T> totalVariance,
    bool doValidate = true
);

template <std::floating_point T>
core::Matrix<T> volFromTotalVariance(
    const std::span<const T> t,
    const core::Matrix<T>& totalVariance,
    bool doValidate = true
);

template <std::floating_point T> T variance(T vol, bool doValidate = true);

template <std::floating_point T>
void variance(std::span<T> out, std::span<const T> vol, bool doValidate = true);

template <std::floating_point T>
core::Matrix<T> variance(const core::VolSurface<T>& volSurface, bool doValidate = true);

template <std::floating_point T>
T atmParameter(
    std::span<const T> parameters,
    std::span<const T> logKF,
    bool doValidate = true
);

template <std::floating_point T> T impliedVol(T callPrice, T t, T r, T q, T S, T K);

namespace detail
{
double impliedVolJackelCall(double callPrice, double t, double dF, double F, double K);
}
} // namespace uv::math::vol

#include "Math/Functions/Detail/Volatility.inl"