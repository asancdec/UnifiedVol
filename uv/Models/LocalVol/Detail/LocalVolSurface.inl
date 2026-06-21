// SPDX-License-Identifier: Apache-2.0

#include "Math/LinearAlgebra/VectorOps.hpp"

namespace uv::models::local_vol
{

template <std::floating_point T>
LocalVolSlice<T>::LocalVolSlice(std::span<const T> cPoints, std::span<const T> knots)
    : interpolator_(exponentiateControlPoints(cPoints), knots)
{
}

template <std::floating_point T>
Vector<T> LocalVolSlice<T>::exponentiateControlPoints(std::span<const T> cPoints) const
{
    return math::linear_algebra::exponential(cPoints);
}

template <std::floating_point T>
void LocalVolSlice<T>::localVarInplace(std::span<T> out, std::span<const T> x) const
{
    interpolator_.evalInplace(out, x);
}

template <std::floating_point T>
Vector<T> LocalVolSlice<T>::localVar(std::span<const T> x) const
{
    return interpolator_.eval(x);
}

template <std::floating_point T>
void LocalVolSlice<T>::localVolInplace(std::span<T> out, std::span<const T> x) const
{
    localVarInplace(out, x);
    math::linear_algebra::squareRootInplace<T>(out, out);
}

template <std::floating_point T>
Vector<T> LocalVolSlice<T>::localVol(std::span<const T> x) const
{
    Vector<T> out{localVar(x)};
    math::linear_algebra::squareRootInplace<T>(out, out);

    return out;
}

} // namespace uv::models::local_vol
