// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Calibrate/Detail/Contexts.hpp"

namespace uv::models::svi::detail
{

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