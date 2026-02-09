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