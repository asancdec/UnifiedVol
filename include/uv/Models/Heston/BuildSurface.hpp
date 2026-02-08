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
#include "Core/VolSurface.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Price/Config.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <concepts>

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N = price::detail::defaultNodes>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const calibrate::Config& config = {}
);

template <std::floating_point T, std::size_t N>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const price::Pricer<T, N>& pricer
);

} // namespace uv::models::heston

#include "Models/Heston/Detail/BuildSurface.inl"
