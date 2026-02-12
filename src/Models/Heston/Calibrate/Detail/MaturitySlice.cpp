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

#include "Models/Heston/Calibrate/Detail/MaturitySlice.hpp"
#include "Base/Macros/Require.hpp"
#include "Math/Functions/Volatility.hpp"

namespace uv::models::heston::calibrate::detail
{

MaturitySlice::MaturitySlice(std::size_t capacity) noexcept
{
    K.reserve(capacity);
    mkt.reserve(capacity);
    w.reserve(capacity);
}

Vector<MaturitySlice> makeSlices(
    std::span<const double> TEST,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& callPrice,
    const opt::cost::WeightATM<double>& weightATM
)
{
    validateInputs(TEST, discountFactors, forwards, strikes, callPrice);

    const std::size_t numStrikes{strikes.size()};
    const std::size_t numMaturities{TEST.size()};

    Vector<MaturitySlice> out;
    out.reserve(numMaturities);

    Vector<double> bufferWeights(numStrikes);
    Vector<double> bufferLogKF(numStrikes);

    for (std::size_t i = 0; i < numMaturities; ++i)
    {
        const double F{forwards[i]};

        out.emplace_back(MaturitySlice{numStrikes});
        MaturitySlice& s = out.back();

        s.t = TEST[i];
        s.dF = discountFactors[i];
        s.F = F;

        std::span<const double> callPriceRow{callPrice[i]};

        math::vol::logKF<double>(bufferLogKF, F, strikes, true);
        opt::cost::weightsATM<double>(bufferLogKF, weightATM, bufferWeights);

        for (std::size_t j = 0; j < numStrikes; ++j)
        {
            s.K.push_back(strikes[j]);
            s.mkt.push_back(callPriceRow[j]);
            s.w.push_back(bufferWeights[j]);
        }
    }

    return out;
}

void validateInputs(
    const std::span<const double> maturities,
    const std::span<const double> discountFactors,
    const std::span<const double> forwards,
    const std::span<const double> strikes,
    const core::Matrix<double>& callPrice
)
{
    UV_REQUIRE_NON_EMPTY(maturities);
    UV_REQUIRE_NON_EMPTY(discountFactors);
    UV_REQUIRE_NON_EMPTY(forwards);
    UV_REQUIRE_NON_EMPTY(strikes);

    UV_REQUIRE_FINITE(maturities);
    UV_REQUIRE_FINITE(discountFactors);
    UV_REQUIRE_FINITE(forwards);
    UV_REQUIRE_FINITE(strikes);

    UV_REQUIRE_POSITIVE(maturities);
    UV_REQUIRE_POSITIVE(discountFactors);

    UV_REQUIRE_SAME_SIZE(maturities, forwards);
    UV_REQUIRE_SAME_SIZE(maturities, discountFactors);
    UV_REQUIRE_SAME_SIZE(maturities, callPrice.rows());
    UV_REQUIRE_SAME_SIZE(strikes, callPrice.cols());

    for (std::size_t i{0}; i < maturities.size(); ++i)
    {
        std::span<const double> callMRow{callPrice[i]};

        UV_REQUIRE_NON_EMPTY(callMRow);
        UV_REQUIRE_FINITE(callMRow);
        UV_REQUIRE_POSITIVE(callMRow);
    }
}
} // namespace uv::models::heston::calibrate::detail