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

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

namespace uv::math::opt::ceres
{

struct Config
{

    unsigned maxEval;
    double functionTol;
    double paramTol;
    double gradientTol;

    double lossScale;

    std::vector<std::string_view> paramNames;

    bool verbose{false};
};

} // namespace uv::math::opt::ceres