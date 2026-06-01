// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Detail/SliceData.hpp"
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
} // namespace uv::models::svi::detail