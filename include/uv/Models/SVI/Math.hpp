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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <array>
#include <concepts>

namespace uv::models::svi
{

template <std::floating_point T>
T totalVariance(T a, T b, T rho, T m, T sigma, T k) noexcept;

template <std::floating_point T> T gk(T a, T b, T rho, T m, T sigma, T k) noexcept;

namespace detail
{
template <std::floating_point T>
T aParam(T atmTotalVariance, T b, T rho, T m, T sigma) noexcept;

} // namespace detail

} // namespace uv::models::svi

#include "Models/SVI/Detail/Math.inl"