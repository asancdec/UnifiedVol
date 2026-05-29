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

#include <cmath>

namespace uv::models::svi
{
template <std::floating_point T>
T totalVariance(T a, T b, T rho, T m, T sigma, T k) noexcept
{
    const T x{k - m};
    const T R{std::sqrt(std::fma(x, x, sigma * sigma))};
    return std::fma(b, (rho * x + R), a);
}

template <std::floating_point T> T gk(T a, T b, T rho, T m, T sigma, T k) noexcept
{
    T x{k - m};
    T sigmaSquared{sigma * sigma};

    T R{std::sqrt(std::fma(x, x, sigmaSquared))};
    T invR{1.0 / R};

    T wk{std::fma(b, (rho * x + R), a)};
    T wkInv{1.0 / wk};

    T wkD1{b * (rho + x * invR)};
    T wkD1Squared{wkD1 * wkD1};

    T invR2{invR * invR};
    T invRCubed{invR2 * invR};

    T wkD2{b * sigmaSquared * invRCubed};

    T A{1.0 - 0.5 * k * wkD1 * wkInv};
    T B{wkInv + 0.25};

    return (A * A) - 0.25 * wkD1Squared * B + wkD2 * 0.5;
}
} // namespace uv::models::svi

namespace uv::models::svi::detail
{
template <std::floating_point T>
T aParam(T atmTotalVariance, T b, T rho, T m, T sigma) noexcept
{
    return atmTotalVariance - b * (-rho * m + std::sqrt(m * m + sigma * sigma));
}

} // namespace uv::models::svi::detail