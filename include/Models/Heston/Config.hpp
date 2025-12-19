// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Config.hpp
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

#include <limits>
#include <concepts>

namespace uv::models::heston
{
    /**
     * @brief Configuration parameters for the Heston Fourier / contour-shift pricer.
     *
     * Stores damping parameters used to ensure integrability of the pricing integral,
     * plus a small numerical epsilon used for stability guards.
     */
    template<std::floating_point T>
    struct Config
    {
        T alphaItm;                                    ///< Damping parameter when ln(F/K) >= 0 (ITM region).
        T alphaOtm;                                    ///< Damping parameter when ln(F/K) < 0  (OTM region).
        T eps{ std::numeric_limits<T>::epsilon() };    ///< Numerical epsilon (defaults to machine epsilon).
    };
}
