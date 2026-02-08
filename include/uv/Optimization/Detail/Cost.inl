// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 �lvaro S�nchez de Carlos
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

#include "Base/Macros/Require.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>

namespace uv::opt::cost
{
template <std::floating_point T>

void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate
)
{
    if (doValidate)
    {

        detail::validateWeightsATM<T>(logKF, params, out);
    }

    const T wATMMinusOne{params.wATM - 1.0};
    const T invk0{1.0 / params.k0};

    for (std::size_t i{0}; i < logKF.size(); ++i)
    {
        const T z{logKF[i] * invk0};

        out[i] = std::sqrt(1.0 + wATMMinusOne * std::exp(-(z * z)));
    }
}
} // namespace uv::opt::cost

namespace uv::opt::cost::detail
{
template <std::floating_point T>
void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
)
{
    const T wATM{params.wATM};
    const T k0{params.k0};

    UV_REQUIRE_NON_EMPTY(out);
    UV_REQUIRE_NON_EMPTY(logKF);

    UV_REQUIRE_FINITE(logKF);
    UV_REQUIRE_FINITE(wATM);
    UV_REQUIRE_FINITE(k0);

    UV_REQUIRE_SAME_SIZE(out, logKF);

    UV_REQUIRE_EQUAL_OR_GREATER(1.0, Watm);
    UV_REQUIRE_NON_NEGATIVE(k0);
}
} // namespace uv::opt::cost::detail