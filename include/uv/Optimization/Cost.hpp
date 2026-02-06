// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

namespace uv::opt::cost
{
template <std::floating_point T> struct WeightATM
{
    T wATM{1.0};
    T k0{};
};

template <std::floating_point T>
void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate = true
);
namespace detail
{

template <std::floating_point T>
void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
);
}
} // namespace uv::opt::cost

#include <Optimization/Detail/Cost.inl>