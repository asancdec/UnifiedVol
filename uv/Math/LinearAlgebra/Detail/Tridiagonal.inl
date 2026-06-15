// SPDX-License-Identifier: Apache-2.0

namespace uv::math::linear_algebra
{
template <std::floating_point T, std::size_t N> void thomasSolve(
    std::span<T, N> x,
    std::span<const T, N> upper,
    std::span<const T, N> middle,
    std::span<const T, N> lower,
    std::span<T, N> scratch
) noexcept
{
    static_assert(N >= 2, "thomasSolve: N must be >= 2");
    constexpr std::size_t last{N - 1};

    const T middleInv{T{1} / middle[0]};

    scratch[0] = upper[0] * middleInv;
    x[0] *= middleInv;

    for (std::size_t i{1}; i < N; ++i)
    {
        const std::size_t previous{i - 1};
        const T lowerI{lower[i]};
        const T xI{x[i]};
        const T invDenom{T{1} / (middle[i] - lowerI * scratch[previous])};

        if (i < last)
        {
            scratch[i] = upper[i] * invDenom;
        }
        x[i] = (xI - lowerI * x[previous]) * invDenom;
    }

    for (std::size_t i{last}; i-- > 0;)
    {
        x[i] -= scratch[i] * x[i + 1];
    }
}

} // namespace uv::math::linear_algebra
