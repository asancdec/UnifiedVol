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

#include <concepts>

namespace uv::models::heston
{
    /**
     * @brief Heston stochastic volatility model parameters.
     *
     * The Heston model dynamics are commonly written as:
     *  - dS_t = (r - q) S_t dt + sqrt(v_t) S_t dW^S_t
     *  - dv_t = kappa (theta - v_t) dt + sigma sqrt(v_t) dW^v_t
     *  - corr(dW^S_t, dW^v_t) = rho
     *
     * This struct stores the five scalar parameters required by the model.
     */
    template <std::floating_point T>
    struct Params
    {
        T kappa;  ///< Mean reversion speed of variance (kappa > 0).
        T theta;  ///< Long-run variance level (theta >= 0).
        T sigma;  ///< Volatility of variance ("vol-of-vol", sigma >= 0).
        T rho;    ///< Correlation between spot and variance Brownian motions (rho in [-1, 1]).
        T v0;     ///< Initial variance at t = 0 (v0 >= 0).
    };
}