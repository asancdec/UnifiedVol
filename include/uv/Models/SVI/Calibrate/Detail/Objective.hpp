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

#include "Base/Alias.hpp"
#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <concepts>
#include <cstddef>

namespace uv::models::svi::detail
{

template <opt::nlopt::Algorithm Algo>
void setMinObjective(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const ObjectiveContexts& ctx
) noexcept;

double objectiveThunk(unsigned /*n*/, const double* x, double* grad, void* data) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Objective.inl"
