// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Bindings.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
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


#include "Math/Interpolation.hpp"
#include "Utils/Types.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(interp, m)
{
    m.def(
        "interpolate_cubic_hermite",
        [](const std::vector<long double> x,
            const std::vector<long double>& xs,
            const std::vector<long double>& ys,
            const std::vector<long double>& dydx)
        {
            return uv::math::interpolateCubicHermiteSpline<long double>(x, xs, ys, dydx);
        },
        R"pbdoc(
        Cubic Hermite interpolation with precomputed slopes.

        Given node locations `xs`, function values `ys`, and first-derivative estimates
        `dydx` at each node, this function evaluates the C¹-continuous cubic Hermite
        interpolant at a collection of query points `x`.

        On each interval [xs[i], xs[i+1]], the interpolant is a cubic polynomial

            H(t) = y[i]
                   + d[i]  * Δx
                   + c2[i] * (Δx)²
                   + c3[i] * (Δx)³

        where Δx = t - xs[i], and the coefficients are computed locally from the node
        spacing and slopes. No global system is solved.

        Extrapolation outside the data range is clamped:

        * For x < xs[0], the value ys[0] is returned.
        * For x > xs[-1], the value ys[-1] is returned.

        Parameters
        ----------
        x : Sequence[float]
            Query points at which to evaluate the interpolant.
        xs : Sequence[float]
            Strictly increasing x-grid (knot locations).
        ys : Sequence[float]
            Function values corresponding to `xs`.
        dydx : Sequence[float]
            First-derivative estimates at each x in `xs`.

        Returns
        -------
        Sequence[float]
            Interpolated values at each query point in `x`.

        Notes
        -----
        - Requires len(xs) == len(ys) == len(dydx) >= 2.
        - Local method: each interval is computed independently; no global solve.
        - Function value and first derivative are continuous across all knots.
        - Complexity: O(N) preprocessing + O(M log N) evaluation for M = len(x).

        References
        ----------
        - Monotone cubic interpolation:
          https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
        )pbdoc"
        );

    m.def(
        "pchip_derivatives",
        [](const std::vector<long double>& xs,
            const std::vector<long double>& ys)
        {
            return uv::math::pchipDerivatives<long double>(xs, ys);
        },
        R"pbdoc(
        Compute PCHIP (monotone cubic Hermite) derivative estimates
        for a 1D sequence of points (xs, ys).

        Given strictly increasing knot positions xs and corresponding function
        values ys, this returns first-derivative estimates dydx suitable for use
        in cubic Hermite interpolation.

        This slope selection reproduces the MATLAB PCHIP / Fritsch–Carlson
        shape-preserving cubic interpolation and guarantees:

        * C¹ continuity
        * No spurious overshoots
        * Local monotonicity preservation
        * Endpoint slope correction and clamping

        Parameters
        ----------
        xs : Sequence[float]
            Strictly increasing x-grid (knot locations).
        ys : Sequence[float]
            Function values at each knot.

        Returns
        -------
        list[float]
            Derivative estimates dydx at each knot.

        Notes
        -----
        - Requires len(xs) == len(ys) >= 2.
        - This does **not** perform interpolation — use `interpolate_cubic_hermite`
          with xs, ys, dydx to evaluate the spline.
        )pbdoc"
    );
}