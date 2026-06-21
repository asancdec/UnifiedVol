// SPDX-License-Identifier: Apache-2.0

#include "Models/Heston/Calibrate/Detail/MaturitySlice.hpp"
#include "Base/Macros/Require.hpp"
#include "Math/Functions/Volatility.hpp"

namespace uv::models::heston::calibrate::detail
{

MaturitySlice::MaturitySlice(std::size_t capacity)
{
    K.reserve(capacity);
    vol.reserve(capacity);
    w.reserve(capacity);
}

Vector<MaturitySlice> makeSlices(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& vol,
    const opt::cost::WeightATM<double>& weightATM
)
{
    validateInputs(maturities, discountFactors, forwards, strikes, vol);

    const std::size_t numStrikes{strikes.size()};
    const std::size_t numMaturities{maturities.size()};

    Vector<MaturitySlice> out;
    out.reserve(numMaturities);

    Vector<double> bufferWeights(numStrikes);
    Vector<double> bufferLogKF(numStrikes);

    for (std::size_t i = 0; i < numMaturities; ++i)
    {
        const double F{forwards[i]};

        out.emplace_back(numStrikes);
        MaturitySlice& s = out.back();

        s.t = maturities[i];
        s.dF = discountFactors[i];
        s.F = F;

        std::span<const double> volRow{vol[i]};

        math::vol::logKF<double>(bufferLogKF, F, strikes, true);
        opt::cost::weightsATM<double>(bufferLogKF, weightATM, bufferWeights);

        for (std::size_t j = 0; j < numStrikes; ++j)
        {
            s.K.push_back(strikes[j]);
            s.vol.push_back(volRow[j]);
            s.w.push_back(bufferWeights[j]);
        }
    }

    return out;
}

void validateInputs(
    const std::span<const double> maturities,
    const std::span<const double> discountFactors,
    const std::span<const double> forwards,
    const std::span<const double> strikes,
    const core::Matrix<double>& vol
)
{
    REQUIRE_NON_EMPTY(maturities);
    REQUIRE_NON_EMPTY(discountFactors);
    REQUIRE_NON_EMPTY(forwards);
    REQUIRE_NON_EMPTY(strikes);

    REQUIRE_FINITE(maturities);
    REQUIRE_FINITE(discountFactors);
    REQUIRE_FINITE(forwards);
    REQUIRE_FINITE(strikes);

    REQUIRE_POSITIVE(maturities);
    REQUIRE_POSITIVE(discountFactors);

    REQUIRE_SAME_SIZE(maturities, forwards);
    REQUIRE_SAME_SIZE(maturities, discountFactors);
    REQUIRE_SAME_SIZE(maturities, vol.rows());
    REQUIRE_SAME_SIZE(strikes, vol.cols());

    for (std::size_t i{0}; i < maturities.size(); ++i)
    {
        std::span<const double> volRow{vol[i]};

        REQUIRE_NON_EMPTY(volRow);
        REQUIRE_FINITE(volRow);
        REQUIRE_POSITIVE(volRow);
    }
}
} // namespace uv::models::heston::calibrate::detail
