// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/DevStatus.hpp"
#include "Base/Macros/Require.hpp"

#include <cmath>

namespace uv::core
{
template <std::floating_point T>
Curve<T>::Curve(T continuouslyCompoundedRate, std::span<const T> maturities)

    : numMaturities_(maturities.size()),
      maturities_(maturities.begin(), maturities.end()),
      discountFactors_(numMaturities_)

{
    UV_REQUIRE_NON_EMPTY(maturities_);
    UV_REQUIRE_FINITE(continuouslyCompoundedRate);
    UV_REQUIRE_FINITE(maturities_);
    UV_REQUIRE_NON_NEGATIVE(maturities_);
    UV_REQUIRE_STRICTLY_INCREASING(maturities_);

    for (std::size_t i{0}; i < numMaturities_; ++i)
    {
        discountFactors_[i] = std::exp(-continuouslyCompoundedRate * maturities_[i]);
    }
}

template <std::floating_point T>
T Curve<T>::interpolateDF(T maturity, bool doValidate) const
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(maturity);
        UV_REQUIRE_NON_NEGATIVE(maturity);
    }

    for (std::size_t i{0}; i < numMaturities_; ++i)
    {
        if (maturity == maturities_[i])
        {
            return discountFactors_[i];
        }
    }

    UV_NOT_IMPLEMENTED("Curve interpolation");
}

template <std::floating_point T>
Vector<T> Curve<T>::interpolateDF(std::span<const T> maturities, bool doValidate) const
{
    if (doValidate)
    {
        UV_REQUIRE_FINITE(maturities);
        UV_REQUIRE_NON_NEGATIVE(maturities);
    }

    const std::size_t n{maturities.size()};

    Vector<T> out;
    out.resize(n);

    for (std::size_t i{0}; i < n; ++i)
    {
        out[i] = interpolateDF(maturities[i], false);
    }

    return out;
}

} // namespace uv::core