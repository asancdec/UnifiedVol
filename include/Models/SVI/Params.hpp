// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Params.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

#include "Utils/Types.hpp"

namespace uv::models::svi
{
    /**
     * @brief SVI parameters for one maturity.
     *
     * Represents a single calibrated SVI slice describing the
     * total variance smile at a given maturity.
     *
     * The parameter `a` is typically fixed from the ATM total variance
     * and not optimized directly.
     */
    struct Params
    {
        Real T;       ///< Maturity (years)
        Real a;       ///< ATM total variance level
        Real b;       ///< Smile amplitude
        Real rho;     ///< Skew parameter (-1 < rho < 1)
        Real m;       ///< Horizontal shift
        Real sigma;   ///< Smile curvature / smoothness
    };
}