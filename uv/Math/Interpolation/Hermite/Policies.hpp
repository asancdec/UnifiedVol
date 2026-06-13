// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <span>

namespace uv::math::interp::hermite
{

template <class D, class T>
concept HasDerivatives = requires(
    const D& d,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> dydx,
    bool doValidate
) {
    { d.derivatives(xs, ys, dydx, doValidate) } -> std::same_as<void>;
};

template <class D, class T>
concept HasEvaluate = requires(
    const D& d,
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
) {
    { d.evaluate(x, xs, ys, dydx, y, doValidate) } -> std::same_as<void>;
};

template <class Derived, std::floating_point T> struct Derivatives
{

    Vector<T>
    operator()(std::span<const T> xs, std::span<const T> ys, bool doValidate = true) const
    requires HasDerivatives<Derived, T>;

    void operator()(
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> dydx,
        bool doValidate = true
    ) const
    requires HasDerivatives<Derived, T>;
};

template <std::floating_point T> struct PchipDerivatives
    : Derivatives<PchipDerivatives<T>, T>
{

    void derivatives(
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> dydx,
        bool doValidate
    ) const;
};

template <class Derived, std::floating_point T> struct Evaluate
{

    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasEvaluate<Derived, T>;

    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate = true
    ) const
    requires HasEvaluate<Derived, T>;

    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasEvaluate<Derived, T>;
};

template <std::floating_point T> struct HermiteEval : Evaluate<HermiteEval<T>, T>
{

    void evaluate(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate
    ) const;
};

namespace detail
{

template <std::floating_point T> void hermiteSplineInterp(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
);

template <std::floating_point T> void pchipDerivatives(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> dydx,
    bool doValidate
);

template <std::floating_point T> T pchipEndpointSlope(T h1, T h2, T S1, T S2) noexcept;

template <std::floating_point T> void validateInputsDerivatives(
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx
);

template <std::floating_point T> void validateInputsEvaluate(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<const T> y
);

} // namespace detail
} // namespace uv::math::interp::hermite

#include "Math/Interpolation/Hermite/Detail/Policies.inl"
