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

#include <Models/SVI/Calibrate/Detail/Contexts.hpp>
#include <Models/SVI/Math.hpp>

namespace uv::models::svi::detail
{

ObjectiveContexts::ObjectiveContexts(
    std::span<const double> logKF,
    std::span<const double> totalVariance,
    double atmTotalVariance
) noexcept
    : k(logKF.data()),
      wM(totalVariance.data()),
      n(logKF.size()),
      atmTotalVariance(atmTotalVariance)
{
}
} // namespace uv::models::svi::detail