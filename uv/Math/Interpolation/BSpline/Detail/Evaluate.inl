// SPDX-License-Identifier: Apache-2.0

#include "Base/Types.hpp"

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::interp::bspline::detail
{

template <std::floating_point T, int K> Vector<T>
eval(std::span<const T> x, std::span<const T> knots, std::span<const T> cPoints) noexcept
{
    Vector<T> out(x.size());

    evalInplace<T, K>(out, x, knots, cPoints);

    return out;
}

template <std::floating_point T, int K> void evalInplace(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::span<const T> cPoints
) noexcept
{
    static_assert(K >= 0, "Spline degree must be non-negative");

    const std::size_t xSize{x.size()};

    std::fill(out.begin(), out.end(), T{0});

    Vector<T> basisValues(xSize);

    for (std::size_t i{0}; i < cPoints.size(); ++i)
    {
        coxDeBoor<T, K>(basisValues, x, knots, i);

        const T cPoint{cPoints[i]};
        for (std::size_t j{0}; j < xSize; ++j)
        {
            out[j] += basisValues[j] * cPoint;
        }
    }

    for (std::size_t j{0}; j < xSize; ++j)
    {
        if (x[j] == knots.back())
        {
            out[j] = cPoints.back();
        }
    }
}

template <std::floating_point T, int K> void coxDeBoor(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::size_t idx
) noexcept
{
    static_assert(K >= 0, "Basis degree must be non-negative");

    if constexpr (K == 0)
    {
        basisZero(out, x, knots, idx);
    }
    else
    {
        const std::size_t xSize{x.size()};

        std::fill(out.begin(), out.end(), T{0});

        const T denomLeft{knots[idx + K] - knots[idx]};
        const T denomRight{knots[idx + K + 1] - knots[idx + 1]};

        if (denomLeft != T{0})
        {
            Vector<T> left(xSize);
            coxDeBoor<T, K - 1>(left, x, knots, idx);

            const T invLeft{T{1} / denomLeft};
            const T knotLower{knots[idx]};

            for (std::size_t j{0}; j < xSize; ++j)
            {
                out[j] += (x[j] - knotLower) * invLeft * left[j];
            }
        }

        if (denomRight != T{0})
        {
            Vector<T> right(xSize);
            coxDeBoor<T, K - 1>(right, x, knots, idx + 1);

            const T invRight{T{1} / denomRight};
            const T knotUpper{knots[idx + K + 1]};

            for (std::size_t j{0}; j < xSize; ++j)
            {
                out[j] += (knotUpper - x[j]) * invRight * right[j];
            }
        }
    }
}

template <std::floating_point T> void basisZero(
    std::span<T> out,
    std::span<const T> x,
    std::span<const T> knots,
    std::size_t idx
) noexcept
{
    const T knotLower{knots[idx]};
    const T knotUpper{knots[idx + 1]};

    for (std::size_t i{0}; i < out.size(); ++i)
    {
        const T elem{x[i]};

        out[i] = (knotLower <= elem && elem < knotUpper) ? T{1} : T{0};
    }
}

} // namespace uv::math::interp::bspline::detail
