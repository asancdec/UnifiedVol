// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Policies.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Interpolation policies and definitions
 *
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
#include "Utils/Aux/Errors.hpp"

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
    // ---------- Validate ----------

    if (doValidate)
    {
        validateInputsEvaluate<T>(x, xs, ys, dydx, y);
    }

    // ---------- Calculate ----------

    const std::size_t xsSize{xs.size()};
    const std::size_t xSize{x.size()};

    const std::size_t numSteps{xsSize - 1};
    Vector<T> h(numSteps); // Step sizes
    Vector<T> S(numSteps); // Secant slopes

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T hi{xs[i + 1] - xs[i]};
        h[i] = hi;
        S[i] = (ys[i + 1] - ys[i]) / hi;
    }

    // Polynomial coefficients

    Vector<T> c2s(numSteps); // Second order
    Vector<T> c3s(numSteps); // Third order

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T c1{dydx[i]}; // Tangent slope
        const T m{S[i]};     // Secant slope
        const T invH{1.0 / h[i]};
        const T common{c1 + dydx[i + 1] - 2.0 * m};

        c2s[i] = (m - c1 - common) * invH;
        c3s[i] = (common * invH * invH);
    }

    // Interpolate
    const T xsMin{xs.front()};
    const T xsMax{xs.back()};

    for (std::size_t i{0}; i < xSize; ++i)
    {
        // Out of bounds
        const T xi{x[i]};

        // Left bounds
        if (xi <= xsMin)
        {
            // Flat extrapolation
            y[i] = ys.front();
            continue;
        }

        // Right bounds
        if (xi >= xsMax)
        {
            // Flat extrapolation
            y[i] = ys.back();
            continue;
        }

        // ---------- Search ----------

        auto it = std::upper_bound(xs.begin(), xs.end(), xi);
        std::size_t idx = static_cast<std::size_t>(it - xs.begin()) - 1;

        // ---------- Evaluate ----------

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
    // ---------- Validate ----------

    if (doValidate)
        validateInputsDerivatives<T>(xs, ys, dydx);

    // ---------- Special cases  ----------

    const std::size_t xsSize{xs.size()};

    if (xsSize == 2)
    {
        // Single secant slope
        const T S{(ys.back() - ys.front()) / (xs.back() - xs.front())};

        dydx.front() = S;
        dydx.back() = S;

        // End function
        return;
    }

    // ---------- Calculate ----------

    // Secant slopes

    const std::size_t numSteps{xsSize - 1};

    Vector<T> h(numSteps); // Step sizes
    Vector<T> S(numSteps); // Secant slopes

    for (std::size_t i{0}; i < numSteps; ++i)
    {
        const T hi{xs[i + 1] - xs[i]};
        h[i] = hi;
        S[i] = (ys[i + 1] - ys[i]) / hi;
    }

    // Inner derivatives

    // Set to zero
    std::fill(dydx.begin(), dydx.end(), 0.0);

    for (std::size_t i{1}; i < numSteps; ++i)
    {
        const T S1{S[i - 1]};
        const T S2{S[i]};
        const T S1timesS2{S1 * S2};

        // Check slope signs
        // NOTE: else, the derivative is already set to zero
        if (S1timesS2 > 0.0)
        {
            // Extract
            const T h1{h[i - 1]};
            const T h2{h[i]};

            // Calculate weight
            const T weight{(h1 + 2.0 * h2) / (3.0 * (h1 + h2))};

            // Calculate tangent slope (derivative)
            dydx[i] = S1timesS2 / (weight * S2 + (1.0 - weight) * S1);
        }
    }

    // Left endpoint
    dydx.front() = pchipEndpointSlope<T>(
        h[0], // h1 = xs[1] - xs[0]
        h[1], // h2 = xs[2] - xs[1]
        S[0], // S1 = (ys[1] - ys[0]) / h1
        S[1]  // S2 = (ys[2] - ys[1]) / h2
    );

    // Right endpoint
    dydx.back() = pchipEndpointSlope<T>(
        h[numSteps - 1], // h1 = xs[n-1] - xs[n-2]
        h[numSteps - 2], // h2 = xs[n-2] - xs[n-3]
        S[numSteps - 1], // S1 = (ys[n-1] - ys[n-2]) / h1
        S[numSteps - 2]  // S2 = (ys[n-2] - ys[n-3]) / h2
    );
}

template <std::floating_point T>
T pchipEndpointSlope(const T h1, const T h2, const T S1, const T S2) noexcept
{
    // ---------- Calculate ----------

    // Derivative
    T d{((2.0 * h1 + h2) * S1 - h1 * S2) / (h1 + h2)};

    // ---------- Check ----------

    // If d points against the local trend, clamp to zero
    if (std::signbit(d) != std::signbit(S1))
    {
        d = 0.0;
    }
    // If adjacent secant slopes disagree AND magnitude too large, clamp
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
    // ---------- Size ----------

    const std::size_t xsSize{xs.size()};
    const std::size_t ysSize{ys.size()};

    // Throw if sizes do not match
    UV_REQUIRE(
        xsSize == ysSize,
        ErrorCode::InvalidArgument,
        "validateInputs: size mismatch — xs vector has " + std::to_string(xsSize) +
            " elements but ys vector has " + std::to_string(ysSize) + " elements"
    );

    // Throw if dydx size does not match xs/ys
    UV_REQUIRE(
        dydx.size() == xsSize,
        ErrorCode::InvalidArgument,
        "validateInputs: size mismatch — dydx vector has " + std::to_string(dydx.size()) +
            " elements but xs/ys have " + std::to_string(xsSize) + " elements"
    );

    // Throw if one point or less
    UV_REQUIRE(
        xsSize > 1,
        ErrorCode::InvalidArgument,
        "validateInputs: xs/ys must contain at least 2 points, but got " +
            std::to_string(xsSize)
    );

    // ---------- NaNs ----------

    for (std::size_t i = 0; i < xsSize; ++i)
    {
        UV_REQUIRE(
            std::isfinite(xs[i]),
            ErrorCode::InvalidArgument,
            "validateInputs: xs[" + std::to_string(i) + "] = " + std::to_string(xs[i]) +
                " is not finite"
        );

        UV_REQUIRE(
            std::isfinite(ys[i]),
            ErrorCode::InvalidArgument,
            "validateInputs: ys[" + std::to_string(i) + "] = " + std::to_string(ys[i]) +
                " is not finite"
        );
    }

    // ---------- Monotonicity ----------

    bool isMonotonic{true};
    std::size_t violIndex{0};

    for (std::size_t i = 1; i < xsSize; ++i)
    {
        if (xs[i] <= xs[i - 1])
        {
            isMonotonic = false;
            violIndex = i;
            break;
        }
    }

    // Throw if not strictly monotonically increasing
    UV_REQUIRE(
        isMonotonic,
        ErrorCode::InvalidArgument,
        "validateInputs: xs vector must be strictly increasing "
        "(violation at index " +
            std::to_string(violIndex) +
            ": xs[i-1] = " + std::to_string(xs[violIndex - 1]) +
            ", xs[i] = " + std::to_string(xs[violIndex]) + ")"
    );
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
    // ---------- Grid ----------

    validateInputsDerivatives<T>(xs, ys, dydx);

    // ---------- Size ----------

    UV_REQUIRE(
        y.size() == x.size(),
        ErrorCode::InvalidArgument,
        "validateInputsHermiteSpline: y has " + std::to_string(y.size()) +
            " elements but x has " + std::to_string(x.size()) + " elements"
    );

    // ---------- NaNs ----------

    // Throw if NaN derivatives exist
    for (std::size_t i = 0; i < dydx.size(); ++i)
    {
        UV_REQUIRE(
            std::isfinite(dydx[i]),
            ErrorCode::InvalidArgument,
            "validateInputsHermiteSpline: dydx[" + std::to_string(i) +
                "] = " + std::to_string(dydx[i]) + " is not finite"
        );
    }

    // Throw if NaN exist in target grid
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        UV_REQUIRE(
            std::isfinite(x[i]),
            ErrorCode::InvalidArgument,
            "validateInputsHermiteSpline: x[" + std::to_string(i) +
                "] = " + std::to_string(x[i]) + " is not finite"
        );
    }
}
} // namespace uv::math::interp::detail