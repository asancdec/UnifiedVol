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

#include "Base/Types.hpp"
#include "Core/Matrix.hpp"
#include "Optimization/Cost.hpp"

#include <cstddef>
#include <span>

namespace uv::models::heston::calibrate::detail
{

struct MaturitySlice
{
    double t;
    double dF;
    double F;
    Vector<double> K;
    Vector<double> mkt;
    Vector<double> w;

    MaturitySlice() = delete;

    explicit MaturitySlice(std::size_t capacity) noexcept;
};

Vector<MaturitySlice> makeSlices(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& callPrice,
    const opt::cost::WeightATM<double>& weightATM
);

void validateInputs(
    const std::span<const double> maturities,
    const std::span<const double> discountFactors,
    const std::span<const double> forwards,
    const std::span<const double> strikes,
    const core::Matrix<double>& callPrice
);

} // namespace uv::models::heston::calibrate::detail
