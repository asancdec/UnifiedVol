// SPDX-License-Identifier: Apache-2.0

namespace uv::models::heston::calibrate::detail
{

template <typename Policy> void setGuessBounds(opt::ceres::Optimizer<Policy>& optimizer)
{
    optimizer.initialize(initGuess(), lowerBounds(), upperBounds());
}

} // namespace uv::models::heston::calibrate::detail