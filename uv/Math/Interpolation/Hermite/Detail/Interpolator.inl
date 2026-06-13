// SPDX-License-Identifier: Apache-2.0

namespace uv::math::interp::hermite
{
template <std::floating_point T, class DerivPolicy, class EvalPolicy>
void Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    std::span<T> y,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    eval(x, xs, ys, dydx, y, doValidate);
}

template <std::floating_point T, class DerivPolicy, class EvalPolicy>
Vector<T> Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    Vector<T> dydx(xs.size());
    deriv(xs, ys, dydx, doValidate);

    Vector<T> y(x.size());

    (*this)(x, xs, ys, dydx, y, doValidate);

    return y;
}

template <std::floating_point T, class DerivPolicy, class EvalPolicy>
void Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<T> y,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    Vector<T> dydx(xs.size());
    deriv(xs, ys, dydx, doValidate);

    (*this)(x, xs, ys, dydx, y, doValidate);
}

template <std::floating_point T, class DerivPolicy, class EvalPolicy>
Vector<T> Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    std::span<const T> x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    Vector<T> y(x.size());

    (*this)(x, xs, ys, dydx, y, doValidate);

    return y;
}

template <std::floating_point T, class DerivPolicy, class EvalPolicy>
T Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    T x,
    std::span<const T> xs,
    std::span<const T> ys,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    Vector<T> dydx(xs.size());
    deriv(xs, ys, dydx, doValidate);

    const std::array<T, 1> xIn{x};
    std::array<T, 1> y{};

    (*this)(xIn, xs, ys, dydx, y, doValidate);

    return y.front();
}

template <std::floating_point T, class DerivPolicy, class EvalPolicy>
T Interpolator<T, DerivPolicy, EvalPolicy>::operator()(
    T x,
    std::span<const T> xs,
    std::span<const T> ys,
    std::span<const T> dydx,
    bool doValidate
) const
requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
{
    const std::array<T, 1> xIn{x};
    std::array<T, 1> y{};

    (*this)(xIn, xs, ys, dydx, y, doValidate);

    return y.front();
}

} // namespace uv::math::interp::hermite
