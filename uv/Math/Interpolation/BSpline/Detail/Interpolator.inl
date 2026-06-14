#include "Base/Macros/Require.hpp"
#include "Math/Interpolation/BSpline/Detail/Evaluate.hpp"

#include <cstddef>

namespace uv::math::interp::bspline
{

template <std::floating_point T, int K, bool doValidate>
BSpline<T, K, doValidate>::BSpline(std::span<const T> cPoints, std::span<const T> knots)
    : cPoints_(cPoints.begin(), cPoints.end()),
      knots_(knots.begin(), knots.end())
{
    if constexpr (doValidate)
    {
        validate();
    }
}

template <std::floating_point T, int K, bool doValidate>
void BSpline<T, K, doValidate>::evalInplace(std::span<T> out, std::span<const T> x) const
{
    if constexpr (doValidate)
    {
        REQUIRE_SAME_SIZE(out, x);
    }

    detail::evalInplace<T, K>(out, x, knots_, cPoints_);
}

template <std::floating_point T, int K, bool doValidate>
Vector<T> BSpline<T, K, doValidate>::eval(std::span<const T> x) const noexcept
{
    return detail::eval<T, K>(x, knots_, cPoints_);
}

template <std::floating_point T, int K, bool doValidate>
void BSpline<T, K, doValidate>::validate() const
{
    static_assert(K >= 0, "BSpline degree must be non-negative");

    constexpr std::size_t degree{static_cast<std::size_t>(K)};

    REQUIRE_NON_EMPTY(cPoints_);
    REQUIRE_NON_EMPTY(knots_);
    REQUIRE_MIN_SIZE(cPoints_, degree + 1);
    REQUIRE_EQUAL(knots_.size(), cPoints_.size() + degree + 1);

    REQUIRE_FINITE(cPoints_);
    REQUIRE_FINITE(knots_);
    REQUIRE_NON_DECREASING(knots_);
    REQUIRE_LESS(knots_[degree], knots_[knots_.size() - degree - 1]);
}

} // namespace uv::math::interp::bspline
