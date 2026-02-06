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

#include <Models/SVI/Detail/Params.inl>