#include "Base/Macros/Require.hpp"
#include "Math/Interpolation/BSpline/Detail/Evaluate.hpp"

#include <utility>

namespace uv::math::interp::bspline
{

template <std::floating_point T, std::size_t K, bool doValidate>
BSpline<T, K, doValidate>::BSpline(std::span<const T> cPoints, std::span<const T> knots)
    : cPoints_(cPoints.begin(), cPoints.end()),
      knots_(knots.begin(), knots.end())
{
    static_assert(K >= 0, "BSpline degree must be non-negative");

    if constexpr (doValidate)
    {
        validate();
    }
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::setControlPoints(std::span<const T> cPoints)
{
    Vector<T> nextCPoints(cPoints.begin(), cPoints.end());

    if constexpr (doValidate)
    {
        validate(nextCPoints, knots_);
    }

    cPoints_ = std::move(nextCPoints);
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::setKnots(std::span<const T> knots)
{
    Vector<T> nextKnots(knots.begin(), knots.end());

    if constexpr (doValidate)
    {
        validate(cPoints_, nextKnots);
    }

    knots_ = std::move(nextKnots);
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::evalInplace(std::span<T> out, std::span<const T> x) const
{
    if constexpr (doValidate)
    {
        REQUIRE_SAME_SIZE(out, x);
    }

    detail::evalInplace<T, K>(out, x, knots_, cPoints_);
}

template <std::floating_point T, std::size_t K, bool doValidate>
Vector<T> BSpline<T, K, doValidate>::eval(std::span<const T> x) const
{
    return detail::eval<T, K>(x, knots_, cPoints_);
}

template <std::floating_point T, std::size_t K, bool doValidate> void
BSpline<T, K, doValidate>::validate(std::span<const T> cPoints, std::span<const T> knots)
{
    REQUIRE_NON_EMPTY(cPoints);
    REQUIRE_NON_EMPTY(knots);
    REQUIRE_MIN_SIZE(cPoints, K + 1);
    REQUIRE_EQUAL(knots.size(), cPoints.size() + K + 1);

    REQUIRE_FINITE(cPoints);
    REQUIRE_FINITE(knots);
    REQUIRE_NON_DECREASING(knots);
    REQUIRE_LESS(knots[K], knots[knots.size() - K - 1]);
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::validate() const
{
    validate(cPoints_, knots_);
}

} // namespace uv::math::interp::bspline
