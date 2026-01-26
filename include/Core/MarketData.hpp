// SPDX-License-Identifier: Apache-2.0
/*
 * File:        MarketData.hpp
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

#include <concepts>

namespace uv::core
{
/**
 * @brief Container for basic market inputs.
 *
 * Holds the minimal set of continuously compounded market parameters
 * required by pricing and calibration routines.
 *
 * All quantities are expressed under the same measure and conventions:
 * - interest rates and yields are continuously compounded
 * - spot price is expressed in domestic currency units
 *
 * This struct is a passive data holder with no invariants enforced.
 */
template <std::floating_point T> struct MarketData
{
    T r; // Continuously compounded risk-free rate
    T q; // Continuously compounded dividend yield
    T S; // Spot price
};
} // namespace uv::core