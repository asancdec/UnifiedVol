// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Math/Interpolation/Hermite/Policies.hpp"

#include <concepts>
#include <span>

namespace uv::math::interp::hermite
{

template <std::floating_point T, class DerivPolicy, class EvalPolicy> struct Interpolator
{

    using value_type = T;
    using derivatives_type = DerivPolicy;
    using evaluator_type = EvalPolicy;

    [[no_unique_address]] DerivPolicy deriv{};

    [[no_unique_address]] EvalPolicy eval{};

    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        std::span<T> y,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    void operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<T> y,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    Vector<T> operator()(
        std::span<const T> x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;

    T operator()(
        T x,
        std::span<const T> xs,
        std::span<const T> ys,
        std::span<const T> dydx,
        bool doValidate = true
    ) const
    requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>;
};

template <std::floating_point T> using PchipInterpolator =
    Interpolator<T, PchipDerivatives<T>, HermiteEval<T>>;

} // namespace uv::math::interp::hermite

#include "Math/Interpolation/Hermite/Detail/Interpolator.inl"
