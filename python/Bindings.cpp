/**
* Bindings.cpp
* Author: Alvaro Sanchez de Carlos
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
        [](double x,
            const uv::Vector<uv::Real>& xs,
            const uv::Vector<uv::Real>& ys,
            const uv::Vector<uv::Real>& dydx)
        {
            return uv::math::interpolateCubicHermiteSpline<uv::Real>(x, xs, ys, dydx);
        },
        R"pbdoc(
        Cubic Hermite interpolation with precomputed slopes.

        Given node locations `xs`, function values `ys`, and first-derivative estimates
        `dydx` at each node, this function evaluates the C¹-continuous cubic Hermite
        interpolant at a query point `x`.

        On each interval [x[i], x[i+1]], the interpolant is a cubic polynomial

            H(x) = y[i]
                   + d[i] * Δx
                   + c2[i] * (Δx)²
                   + c3[i] * (Δx)³

        where Δx = x - x[i], and the coefficients are computed locally from the node
        spacing and slopes. No global system is solved.

        The interpolation is clamped outside the data range:

        * For x < xs[0], the value ys[0] is returned.
        * For x > xs[-1], the value ys[-1] is returned.

        Parameters
        ----------
        x : float
            Query point at which to evaluate the interpolant.
        xs : Sequence[float]
            Strictly increasing x-grid (knot locations).
        ys : Sequence[float]
            Function values corresponding to `xs`.
        dydx : Sequence[float]
            First-derivative estimates at each x in `xs`.

        Returns
        -------
        float
            Interpolated value at x.

        Notes
        -----
        - Requires len(xs) == len(ys) == len(dydx) >= 2.
        - The method is local: each interval is computed independently.
        - Function value and first derivative are continuous across all knots.

        References
        ----------
        - Monotone cubic interpolation:
          https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
        )pbdoc"
        );
}