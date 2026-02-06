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

#include "Core/Functions.hpp"
#include "Math/PDE/Functions.hpp"
#include "Utils/Aux/Errors.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>

namespace uv::models::localvol
{
template <std::floating_point T, std::size_t NT, std::size_t NX, class Interpolator>
template <typename F>
Pricer<T, NT, NX, Interpolator>::Pricer(F&& payoff, T xBound)
    : xGrid_(core::generateGrid<T, NX>(-xBound, xBound)),
      cInit_(math::pde::andreasenHugeInit<T, NX>(xGrid_, std::forward<F>(payoff))),
      xGridInner_(xGrid_.data() + 1, NX - 2),
      cInner_(c_.data() + 1, NX - 2)
{

    UV_REQUIRE(
        xBound > 0.0,
        ErrorCode::InvalidArgument,
        "andreasenHugeInit: xBound must be positive"
    );

    const T dX{(2.0 * xBound) / (NX - 1)};
    const T invDX{1.0 / dX};
    const T invDXDivTwo{0.5 * invDX};
    const T invDXSquared{invDX * invDX};
    const T invDXSquaredTimesTwo{invDXSquared * 2.0};
    const T lowerFactor{invDXSquared + invDXDivTwo};
    const T upperFactor{invDXSquared - invDXDivTwo};

    const T rL{std::exp(-dX)};
    const T rR{std::exp(dX)};

    const T aL{1.0 + rL};
    const T bL{-rL};
    const T aR{1.0 + rR};
    const T bR{-rR};

    ahCache_.invDXSquared = invDXSquared;
    ahCache_.invDXSquaredTimesTwo = invDXSquaredTimesTwo;

    ahCache_.lowerFactor = lowerFactor;
    ahCache_.upperFactor = upperFactor;

    ahCache_.aL = aL;
    ahCache_.bL = bL;
    ahCache_.aR = aR;
    ahCache_.bR = bR;
}

template <std::floating_point T, std::size_t NT, std::size_t NX, class Interpolator>
Vector<T> Pricer<T, NT, NX, Interpolator>::priceNormalized(
    T tenor,
    std::span<const T> logKF,
    const VarianceView<T>& localVarView
)
{

    T dt{tenor / NT};
    T dtDivTwo{0.5 * dt};

    ahCache_.zLower = -dtDivTwo * ahCache_.lowerFactor;
    ahCache_.zMiddle = dtDivTwo * ahCache_.invDXSquaredTimesTwo;
    ahCache_.zUpper = -dtDivTwo * ahCache_.upperFactor;

    c_ = cInit_;

    interpolator_(
        xGridInner_,
        localVarView.logKF,
        localVarView.localVar,
        localVarView.dydx,
        ahCache_.localVar
    );

    math::pde::andreasenHugeSolve<T, NT, NX>(c_, cInner_, ahCache_);

    return interpolator_(logKF, xGrid_, c_);
}

} // namespace uv::models::localvol