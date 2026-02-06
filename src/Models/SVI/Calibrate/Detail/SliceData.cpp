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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#include "Models/SVI/Calibrate/Detail/SliceData.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Math/LinearAlgebra/VectorOps.hpp"

namespace uv::models::svi::detail
{
SliceData::SliceData(std::span<const double> logKF, std::span<const double> totalVariance)
    : atmTotalVariance(math::vol::atmParameter(totalVariance, logKF)),
      logKFMin(math::linear_algebra::minValue(logKF)),
      logKFMax(math::linear_algebra::maxValue(logKF))
{
}
} // namespace uv::models::svi::detail