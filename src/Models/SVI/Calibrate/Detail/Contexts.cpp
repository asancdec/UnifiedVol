// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Math/LinearAlgebra/VectorOps.hpp"

namespace uv::models::svi::detail
{

SliceData::SliceData(std::span<const double> logKF, std::span<const double> totalVariance)
    : atmTotalVariance(math::vol::atmParameter(totalVariance, logKF)),
      logKFMin(math::linear_algebra::minValue(logKF)),
      logKFMax(math::linear_algebra::maxValue(logKF))
{
}

ObjectiveContexts::ObjectiveContexts(
    std::span<const double> logKF,
    std::span<const double> totalVariance,
    double atmTotalVariance
) noexcept
    : k(logKF.data()),
      wM(totalVariance.data()),
      n(logKF.size()),
      atmTotalVariance(atmTotalVariance)
{
}
} // namespace uv::models::svi::detail
