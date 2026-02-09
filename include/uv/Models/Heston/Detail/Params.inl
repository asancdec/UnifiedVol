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

#include "Base/Macros/Require.hpp"

namespace uv::models::heston
{

template <std::floating_point T>
constexpr Params<T>::Params(T kappa_, T theta_, T sigma_, T rho_, T v0_) noexcept
    : kappa{kappa_},
      theta{theta_},
      sigma{sigma_},
      rho{rho_},
      v0{v0_}
{
}

template <std::floating_point T> Params<T>::Params(std::span<const double> p)
{
    UV_REQUIRE_NON_EMPTY(p);
    UV_REQUIRE_FINITE(p);
    UV_REQUIRE_SAME_SIZE(p, std::size_t{5});

    kappa = static_cast<T>(p[0]);
    theta = static_cast<T>(p[1]);
    sigma = static_cast<T>(p[2]);
    rho = static_cast<T>(p[3]);
    v0 = static_cast<T>(p[4]);
}

template <std::floating_point T>
template <std::floating_point U>
constexpr Params<U> Params<T>::as() const noexcept
{
    return Params<U>{
        static_cast<U>(kappa),
        static_cast<U>(theta),
        static_cast<U>(sigma),
        static_cast<U>(rho),
        static_cast<U>(v0)
    };
}

} // namespace uv::models::heston