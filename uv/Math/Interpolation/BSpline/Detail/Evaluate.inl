// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, std::size_t K>
std::size_t findSpan(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept
{
    const std::size_t n{cPoints.size() - 1};

    if (!(x < knots[n + 1]) && !(knots[n + 1] < x))
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
    std::array<T, K + 1> d;

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

            if (denom > T{0} || denom < T{0})
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
T evalOne(T x, std::span<const T> knots, std::span<const T> cPoints) noexcept
{
    const std::size_t n{cPoints.size() - 1};

    const T leftDomain{knots[K]};
    const T rightDomain{knots[n + 1]};

    if (x < leftDomain || x > rightDomain)
    {
        return T{0};
    }

    if (!(x < rightDomain))
    {
        return cPoints.back();
    }

    return evalOneInSpan<T, K>(x, findSpan<T, K>(x, knots, cPoints), knots, cPoints);
}

template <std::floating_point T, std::size_t K> void evalInplace(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept
{
    const std::size_t n{cPoints.size() - 1};
    const T leftDomain{knots[K]};
    const T rightDomain{knots[n + 1]};

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

        if (!(xi < rightDomain))
        {
            out[i] = cPoints.back();
            previous = xi;
            continue;
        }

        if (!haveSpan || xi < previous)
        {
            span = findSpan<T, K>(xi, knots, cPoints);
            haveSpan = true;
        }
        else
        {
            while (span < n && xi >= knots[span + 1])
            {
                ++span;
            }
        }

        out[i] = evalOneInSpan<T, K>(xi, span, knots, cPoints);
        previous = xi;
    }
}

template <std::floating_point T, std::size_t K>
Vector<T> eval(std::span<const T> x, std::span<const T> knots, std::span<const T> cPoints)
{
    Vector<T> out(x.size());

    evalInplace<T, K>(out, x, knots, cPoints);

    return out;
}

} // namespace uv::math::interp::bspline::detail
