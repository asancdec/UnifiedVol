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

#include <Base/Macros/Unreachable.hpp>
#include <Optimization/NLopt/Algorithm.hpp>

#include <nlopt.hpp>

namespace uv::opt::nlopt::detail
{
constexpr ::nlopt::algorithm toNlopt(Algorithm a) noexcept
{
    switch (a)
    {
    case Algorithm::LD_MMA:
        return ::nlopt::LD_MMA;
    case Algorithm::LD_SLSQP:
        return ::nlopt::LD_SLSQP;
    case Algorithm::LN_BOBYQA:
        return ::nlopt::LN_BOBYQA;
    }

    UV_UNREACHABLE("Unknown uv::opt::nlopt::Algorithm value");
}
} // namespace uv::opt::nlopt::detail