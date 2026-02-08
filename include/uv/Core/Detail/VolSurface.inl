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

#include "Base/Macros/Require.hpp"

namespace uv::core
{

template <std::floating_point T>
VolSurface<T>::VolSurface(
    std::span<const T> maturities,
    std::span<const T> forwards,
    std::span<const T> strikes,
    std::span<const T> moneyness,
    const Matrix<T>& vol
)
    : maturities_(maturities.begin(), maturities.end()),
      numMaturities_(maturities_.size()),
      strikes_(strikes.begin(), strikes.end()),
      numStrikes_(strikes_.size()),
      forwards_(forwards.begin(), forwards.end()),
      moneyness_(moneyness.begin(), moneyness.end()),
      vol_(vol)
{
    UV_REQUIRE_NON_EMPTY(maturities_);
    UV_REQUIRE_NON_EMPTY(strikes_);
    UV_REQUIRE_NON_EMPTY(forwards_);
    UV_REQUIRE_NON_EMPTY(moneyness_);

    UV_REQUIRE_FINITE(maturities_);
    UV_REQUIRE_FINITE(strikes_);
    UV_REQUIRE_FINITE(forwards_);

    UV_REQUIRE_NON_NEGATIVE(maturities_);
    UV_REQUIRE_NON_NEGATIVE(moneyness_);

    UV_REQUIRE_STRICTLY_INCREASING(maturities_);
    UV_REQUIRE_STRICTLY_INCREASING(strikes_);
    UV_REQUIRE_STRICTLY_INCREASING(moneyness_);

    UV_REQUIRE_SAME_SIZE(maturities_, forwards_);
    UV_REQUIRE_SAME_SIZE(maturities_, vol_.rows());
    UV_REQUIRE_SAME_SIZE(strikes_, moneyness_);
    UV_REQUIRE_SAME_SIZE(strikes_, vol_.cols());

    for (std::size_t i{0}; i < numMaturities_; ++i)
    {
        std::span<const T> volSlice{vol_[i]};

        UV_REQUIRE_NON_EMPTY(volSlice);
        UV_REQUIRE_FINITE(volSlice);
        UV_REQUIRE_NON_NEGATIVE(volSlice);
    }
}

template <std::floating_point T> std::size_t VolSurface<T>::numMaturities() const noexcept
{
    return numMaturities_;
}

template <std::floating_point T> std::size_t VolSurface<T>::numStrikes() const noexcept
{
    return numStrikes_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::maturities() const noexcept
{
    return maturities_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::forwards() const noexcept
{
    return forwards_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::strikes() const noexcept
{
    return strikes_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::moneyness() const noexcept
{
    return moneyness_;
}

template <std::floating_point T> const Matrix<T>& VolSurface<T>::vol() const noexcept
{
    return vol_;
}

} // namespace uv::core
