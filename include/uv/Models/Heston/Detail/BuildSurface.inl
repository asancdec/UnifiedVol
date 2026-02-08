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
#include "Core/Generate.hpp"
#include "Core/Matrix.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/Heston/Calibrate/CeresAdapter.hpp"

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const calibrate::Config& config
)
{
    models::heston::price::Pricer<T> hestonPricer{};

    auto hestonOptimizer{calibrate::detail::makeOptimizer(config)};

    Params<T> hestonParams{calibrate::calibrate(
        volSurface,
        curve,
        hestonOptimizer,
        opt::cost::WeightATM{.wATM = 8.0, .k0 = 0.3},
        hestonPricer
    )};

    hestonPricer.setParams(hestonParams);

    return buildSurface<T>(volSurface, curve, hestonPricer);
}

template <std::floating_point T, std::size_t N>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const price::Pricer<T, N>& pricer
)
{
    return core::generateVolSurface<T>(
        volSurface,
        math::vol::impliedVol<T>(pricer.callPrice(volSurface, curve), volSurface, curve)
    );
}
} // namespace uv::models::heston
