#include "Base/Macros/Require.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, std::size_t K>
std::size_t findSpan(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept;

template <std::floating_point T, std::size_t K> T evalOneInSpan(
    T x,
    std::size_t span,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept;

template <std::floating_point T, std::size_t K>
void validate(std::span<const T> cPoints, std::span<const T> knots);

} // namespace uv::math::interp::bspline::detail

namespace uv::math::interp::bspline
{

template <std::floating_point T, std::size_t K, bool doValidate>
BSpline<T, K, doValidate>::BSpline(std::span<const T> cPoints, std::span<const T> knots)
    : cPoints_(cPoints.begin(), cPoints.end()),
      knots_(knots.begin(), knots.end())
{
    if constexpr (doValidate)
    {
        detail::validate<T, K>(cPoints_, knots_);
    }
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::setControlPoints(std::span<const T> cPoints)
{
    Vector<T> nextCPoints(cPoints.begin(), cPoints.end());

    if constexpr (doValidate)
    {
        detail::validate<T, K>(nextCPoints, knots_);
    }

    cPoints_ = std::move(nextCPoints);
}

template <std::floating_point T, std::size_t K, bool doValidate>
void BSpline<T, K, doValidate>::setKnots(std::span<const T> knots)
{
    Vector<T> nextKnots(knots.begin(), knots.end());

    if constexpr (doValidate)
    {
        detail::validate<T, K>(cPoints_, nextKnots);
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

    const std::size_t n{cPoints_.size() - 1};
    const T leftDomain{knots_[K]};
    const T rightDomain{knots_[n + 1]};

    bool haveSpan{false};
    T previous{};
    std::size_t span{K};

    for (std::size_t i{0}; i < x.size(); ++i)
    {
        const T xi{x[i]};

        if (xi < leftDomain || xi > rightDomain)
        {
            out[i] = T{0};
            haveSpan = false;
            previous = xi;
            continue;
        }

        // codeql-suppress[cpp/equality-on-floats]
        if (xi == rightDomain)
        {
            out[i] = cPoints_.back();
            previous = xi;
            continue;
        }

        if (!haveSpan || xi < previous)
        {
            span = detail::findSpan<T, K>(xi, knots_, cPoints_);
            haveSpan = true;
        }
        else
        {
            while (span < n && xi >= knots_[span + 1])
            {
                ++span;
            }
        }

        out[i] = detail::evalOneInSpan<T, K>(xi, span, knots_, cPoints_);
        previous = xi;
    }
}

template <std::floating_point T, std::size_t K, bool doValidate>
Vector<T> BSpline<T, K, doValidate>::eval(std::span<const T> x) const
{
    Vector<T> out(x.size());

    evalInplace(out, x);

    return out;
}

} // namespace uv::math::interp::bspline

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, std::size_t K>
std::size_t findSpan(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept
{
    const std::size_t n{cPoints.size() - 1};

    // codeql-suppress[cpp/equality-on-floats]
    if (x == knots[n + 1])
    {
        return n;
    }

    const auto first{knots.begin() + static_cast<std::ptrdiff_t>(K)};
    const auto last{knots.begin() + static_cast<std::ptrdiff_t>(n + 2)};

    return static_cast<std::size_t>(std::upper_bound(first, last, x) - knots.begin() - 1);
}

template <std::floating_point T, std::size_t K> T evalOneInSpan(
    T x,
    const std::size_t span,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept
{
    std::array<T, K + 1> d{};

    for (std::size_t j{0}; j <= K; ++j)
    {
        d[j] = cPoints[span - K + j];
    }

    for (std::size_t r{1}; r <= K; ++r)
    {
        for (std::size_t j{K + 1}; j-- > r;)
        {
            const std::size_t i{span - K + j};

            const T denom{knots[i + K + 1 - r] - knots[i]};

            // codeql-suppress[cpp/equality-on-floats]
            if (denom != T{0})
            {
                const T alpha{(x - knots[i]) / denom};
                d[j] = (T{1} - alpha) * d[j - 1] + alpha * d[j];
            }
            else
            {
                d[j] = d[j - 1];
            }
        }
    }

    return d[K];
}

template <std::floating_point T, std::size_t K>
void validate(std::span<const T> cPoints, std::span<const T> knots)
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

} // namespace uv::math::interp::bspline::detail
