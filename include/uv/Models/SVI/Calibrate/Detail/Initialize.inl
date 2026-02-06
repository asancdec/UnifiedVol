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

namespace uv::models::svi::detail
{
template <std::floating_point T, opt::nlopt::Algorithm Algo>
void setGuessBounds(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const Params<T>* prevParams,
    const SliceData& sliceData
) noexcept
{
    const std::array<double, 4> lowerBounds{detail::lowerBounds(sliceData.logKFMin)};
    const std::array<double, 4> upperBounds{detail::upperBounds(sliceData.logKFMax)};

    if (!prevParams)
    {

        optimizer.setGuessBounds(coldGuess(), lowerBounds, upperBounds);

        return;
    }

    const Params<double> prevParamsD{prevParams->template as<double>()};
    optimizer.setGuessBounds(warmGuess(prevParamsD), lowerBounds, upperBounds);
}
} // namespace uv::models::svi::detail