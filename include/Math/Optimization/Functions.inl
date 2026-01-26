// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-01-26
 *
 * Description:
 *   Optimization helper utilities.
 *
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

#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>

namespace uv::math::opt
{
template <std::floating_point T>
void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate
)
{
    // ---------- Validate ----------

    if (doValidate)
    {

        detail::validateWeightsATM<T>(logKF, params, out);
    }

    // ---------- Precompute ----------

    const T wATMMinusOne{params.wATM - 1.0};
    const T invk0{1.0 / params.k0};

    // ---------- Calculate ----------

    for (std::size_t i{0}; i < logKF.size(); ++i)
    {
        const T z{logKF[i] * invk0};

        out[i] = std::sqrt(1.0 + wATMMinusOne * std::exp(-(z * z)));
    }
}
} // namespace uv::math::opt

namespace uv::math::opt::detail
{
template <std::floating_point T>
void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
)
{
    // ---------- Extract ----------

    std::size_t numStrikes{logKF.size()};

    const T wATM{params.wATM};
    const T k0{params.k0};

    // ---------- Size ----------

    UV_REQUIRE(
        out.size() == numStrikes,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: size mismatch - out must have the same "
        "size as strikes"
    );

    // ---------- Values ----------

    UV_REQUIRE(
        std::isfinite(wATM) && wATM >= 1.0,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: wATM must be finite and >= 1.0"
    );

    UV_REQUIRE(
        std::isfinite(k0) && k0 > 0.0,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: k0 must be finite and > 0"
    );

    for (std::size_t i{0}; i < numStrikes; ++i)
    {
        UV_REQUIRE(
            std::isfinite(logKF[i]),
            ErrorCode::InvalidArgument,
            "validateWeightsATM: logKF must be finite"
        );
    }
}
} // namespace uv::math::opt::detail