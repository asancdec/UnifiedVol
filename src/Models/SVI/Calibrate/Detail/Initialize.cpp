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

#include <Models/SVI/Calibrate/Detail/Initialize.hpp>

namespace uv::models::svi::detail
{
std::array<double, 4> coldGuess() noexcept
{
    return {0.1, -0.5, 0.1, 0.1};
}

std::array<double, 4> warmGuess(const Params<double>& params) noexcept
{
    return {

        params.b,
        params.rho,
        params.m,
        params.sigma
    };
}

std::array<double, 4> lowerBounds(double logKFMin) noexcept
{
    return {0.001, -0.9999, 10.0 * logKFMin, 0.01};
}

std::array<double, 4> upperBounds(double logKFMax) noexcept
{
    return {2.0, 0.9999, 10.0 * logKFMax, 10.0};
}

} // namespace uv::models::svi::detail