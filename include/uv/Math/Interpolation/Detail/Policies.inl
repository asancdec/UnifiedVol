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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>

namespace uv::math::interp
{
template <class Derived, std::floating_point T>
Vector<T> Derivatives<Derived, T>::operator()(
    std::span<const T> xs,
    std::span<const T> ys,
    bool doValidate
) const
requires HasDerivatives<Derived, T>
{
    Vector<T> dydx(xs.size());

    (*this)(xs, ys, dydx, doValidate);

    return dydx;
}

template <class Derived, std::floating_point T>
void Derivatives<Derived, T>::operator()(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> dydx,
    bool doValidate
) const
requires HasDerivatives<Derived, T>
{
    static_cast<const Derived&>(*this).derivatives(xs, ys, dydx, doValidate);
}

template <std::floating_point T>
void PchipDerivatives<T>::derivatives(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> dydx,
    bool doValidate
) const
{
    detail::pchipDerivatives<T>(xs, ys, dydx, doValidate);
}

template <class Derived, std::floating_point T>
Vector<T> Evaluate<Derived, T>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    bool doValidate
) const
requires HasEvaluate<Derived, T>
{
    Vector<T> y(x.size());
    (*this)(x, xs, ys, dydx, y, doValidate);
    return y;
}

template <class Derived, std::floating_point T>
void Evaluate<Derived, T>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
) const
requires HasEvaluate<Derived, T>
{
    static_cast<const Derived&>(*this).evaluate(x, xs, ys, dydx, y, doValidate);
}

template <class Derived, std::floating_point T>
T Evaluate<Derived, T>::operator()(
    T x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    bool doValidate
) const
requires HasEvaluate<Derived, T>
{
    const std::array<T, 1> xIn{x};
    std::array<T, 1> y{};
    (*this)(xIn, xs, ys, dydx, y, doValidate);

    return y.front();
}

template <std::floating_point T>
void HermiteEval<T>::evaluate(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
) const
{
    detail::hermiteSplineInterp<T>(x, xs, ys, dydx, y, doValidate);
}
} // namespace uv::math::interp

namespace uv::math::interp::detail
{
template <std::floating_point T>
void hermiteSplineInterp(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
)
{

    if (doValidate)
    {
        validateInputsEvaluate<T>(x, xs, ys, dydx, y);
    }

    const std::size_t xsSize{xs.size()};
    const std::size_t xSize{x.size()};

    const std::size_t numSteps{xsSize - 1};
    Vector<T> h(numSteps);
    Vector<T> S(numSteps);

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T hi{xs[i + 1] - xs[i]};
        h[i] = hi;
        S[i] = (ys[i + 1] - ys[i]) / hi;
    }

    Vector<T> c2s(numSteps);
    Vector<T> c3s(numSteps);

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T c1{dydx[i]};
        const T m{S[i]};
        const T invH{1.0 / h[i]};
        const T common{c1 + dydx[i + 1] - 2.0 * m};

        c2s[i] = (m - c1 - common) * invH;
        c3s[i] = (common * invH * invH);
    }

    const T xsMin{xs.front()};
    const T xsMax{xs.back()};

    for (std::size_t i{0}; i < xSize; ++i)
    {

        const T xi{x[i]};

        if (xi <= xsMin)
        {

            y[i] = ys.front();
            continue;
        }

        if (xi >= xsMax)
        {

            y[i] = ys.back();
            continue;
        }

        auto it = std::upper_bound(xs.begin(), xs.end(), xi);
        std::size_t idx = static_cast<std::size_t>(it - xs.begin()) - 1;

        const T dx{xi - xs[idx]};

        y[i] = ys[idx] + dydx[idx] * dx + c2s[idx] * dx * dx + c3s[idx] * dx * dx * dx;
    }
}

template <std::floating_point T>
void pchipDerivatives(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> dydx,
    bool doValidate
)
{
    if (doValidate)
        validateInputsDerivatives<T>(xs, ys, dydx);

    const std::size_t xsSize{xs.size()};

    if (xsSize == 2)
    {

        const T S{(ys.back() - ys.front()) / (xs.back() - xs.front())};

        dydx.front() = S;
        dydx.back() = S;

        return;
    }

    const std::size_t numSteps{xsSize - 1};

    Vector<T> h(numSteps);
    Vector<T> S(numSteps);

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T hi{xs[i + 1] - xs[i]};
        h[i] = hi;
        S[i] = (ys[i + 1] - ys[i]) / hi;
    }

    std::fill(dydx.begin(), dydx.end(), 0.0);

    for (std::size_t i{1}; i < numSteps; ++i)
    {
        const T S1{S[i - 1]};
        const T S2{S[i]};
        const T S1timesS2{S1 * S2};

        if (S1timesS2 > 0.0)
        {

            const T h1{h[i - 1]};
            const T h2{h[i]};

            const T weight{(h1 + 2.0 * h2) / (3.0 * (h1 + h2))};

            dydx[i] = S1timesS2 / (weight * S2 + (1.0 - weight) * S1);
        }
    }

    dydx.front() = pchipEndpointSlope<T>(h[0], h[1], S[0], S[1]);

    dydx.back() = pchipEndpointSlope<T>(
        h[numSteps - 1],
        h[numSteps - 2],
        S[numSteps - 1],
        S[numSteps - 2]
    );
}

template <std::floating_point T>
T pchipEndpointSlope(const T h1, const T h2, const T S1, const T S2) noexcept
{

    T d{((2.0 * h1 + h2) * S1 - h1 * S2) / (h1 + h2)};

    if (std::signbit(d) != std::signbit(S1))
    {
        d = 0.0;
    }

    else if ((std::signbit(S1) != std::signbit(S2)) && std::abs(d) > 3.0 * std::abs(S1))
    {
        d = 3.0 * S1;
    }
    return d;
}

template <std::floating_point T>
void validateInputsDerivatives(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx
)
{
    UV_REQUIRE_NON_EMPTY(xs);
    UV_REQUIRE_NON_EMPTY(ys);
    UV_REQUIRE_NON_EMPTY(dydx);

    UV_REQUIRE_SAME_SIZE(xs, ys);
    UV_REQUIRE_SAME_SIZE(xs, dydx);
    UV_REQUIRE_MIN_SIZE(xs, 2);

    UV_REQUIRE_FINITE(xs);
    UV_REQUIRE_FINITE(ys);
    UV_REQUIRE_FINITE(dydx);

    UV_REQUIRE_STRICTLY_INCREASING(xs);
}

template <std::floating_point T>
void validateInputsEvaluate(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<const T> y

)
{
    validateInputsDerivatives<T>(xs, ys, dydx);

    UV_REQUIRE_FINITE(x);
    UV_REQUIRE_NON_EMPTY(x);
    UV_REQUIRE_SAME_SIZE(y, x);
}
} // namespace uv::math::interp::detail