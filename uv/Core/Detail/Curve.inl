// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/DevStatus.hpp"
#include "Base/Macros/Require.hpp"
#include "Math/LinearAlgebra/VectorOps.hpp"

#include <cmath>

namespace uv::core
{
template <std::floating_point T>
Curve<T>::Curve(T continuouslyCompoundedRate, std::span<const T> maturities)

    : numMaturities_(maturities.size()),
      maturities_(maturities.begin(), maturities.end()),
      discountFactors_(numMaturities_)

{
    REQUIRE_NON_EMPTY(maturities_);
    REQUIRE_FINITE(continuouslyCompoundedRate);
    REQUIRE_FINITE(maturities_);
    REQUIRE_NON_NEGATIVE(maturities_);
    REQUIRE_STRICTLY_INCREASING(maturities_);

    discountFactors_ = math::linear_algebra::exponential<T>(
        math::linear_algebra::multiply<T>(maturities_, -continuouslyCompoundedRate)
    );
}

template <std::floating_point T>
T Curve<T>::interpolateDF(T maturity, bool doValidate) const
{
    if (doValidate)
    {
        REQUIRE_FINITE(maturity);
        REQUIRE_NON_NEGATIVE(maturity);
    }

    for (std::size_t i{0}; i < numMaturities_; ++i)
    {
        // codeql-suppress[cpp/equality-on-floats]: exact stored maturity lookup.
        if (maturity == maturities_[i])
        {
            return discountFactors_[i];
        }
    }

    NOT_IMPLEMENTED("Curve interpolation");
}

template <std::floating_point T>
Vector<T> Curve<T>::interpolateDF(std::span<const T> maturities, bool doValidate) const
{
    if (doValidate)
    {
        REQUIRE_FINITE(maturities);
        REQUIRE_NON_NEGATIVE(maturities);
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
