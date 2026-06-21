// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"

namespace uv::core
{

template <std::floating_point T> VolSurface<T>::VolSurface(
    std::span<const T> maturities,
    std::span<const T> forwards,
    std::span<const T> strikes,
    std::span<const T> moneyness,
    const Matrix<T>& vol
)
    : maturities_(maturities.begin(), maturities.end()),
      numMaturities_(maturities_.size()),
      strikes_(strikes.begin(), strikes.end()),
      numStrikes_(strikes_.size()),
      forwards_(forwards.begin(), forwards.end()),
      moneyness_(moneyness.begin(), moneyness.end()),
      vol_(vol)
{
    REQUIRE_NON_EMPTY(maturities_);
    REQUIRE_NON_EMPTY(strikes_);
    REQUIRE_NON_EMPTY(forwards_);
    REQUIRE_NON_EMPTY(moneyness_);

    REQUIRE_FINITE(maturities_);
    REQUIRE_FINITE(strikes_);
    REQUIRE_FINITE(forwards_);

    REQUIRE_NON_NEGATIVE(maturities_);
    REQUIRE_POSITIVE(strikes_);
    REQUIRE_POSITIVE(forwards_);
    REQUIRE_POSITIVE(moneyness_);

    REQUIRE_STRICTLY_INCREASING(maturities_);
    REQUIRE_STRICTLY_INCREASING(strikes_);
    REQUIRE_STRICTLY_INCREASING(moneyness_);

    REQUIRE_SAME_SIZE(maturities_, forwards_);
    REQUIRE_SAME_SIZE(maturities_, vol_.rows());
    REQUIRE_SAME_SIZE(strikes_, moneyness_);
    REQUIRE_SAME_SIZE(strikes_, vol_.cols());

    for (std::size_t i{0}; i < numMaturities_; ++i)
    {
        std::span<const T> volSlice{vol_[i]};

        REQUIRE_NON_EMPTY(volSlice);
        REQUIRE_FINITE(volSlice);
        REQUIRE_NON_NEGATIVE(volSlice);
    }
}

template <std::floating_point T> std::size_t VolSurface<T>::numMaturities() const noexcept
{
    return numMaturities_;
}

template <std::floating_point T> std::size_t VolSurface<T>::numStrikes() const noexcept
{
    return numStrikes_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::maturities() const noexcept
{
    return maturities_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::forwards() const noexcept
{
    return forwards_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::strikes() const noexcept
{
    return strikes_;
}

template <std::floating_point T>
std::span<const T> VolSurface<T>::moneyness() const noexcept
{
    return moneyness_;
}

template <std::floating_point T> const Matrix<T>& VolSurface<T>::vol() const noexcept
{
    return vol_;
}

} // namespace uv::core
