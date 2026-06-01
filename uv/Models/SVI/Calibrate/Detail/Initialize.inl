// SPDX-License-Identifier: Apache-2.0

namespace uv::models::svi::detail
{
template <std::floating_point T, opt::nlopt::Algorithm Algo> void setGuessBounds(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const Params<T>* prevParams,
    const SliceData& sliceData
) noexcept
{
    const std::array<double, 4> lowerBounds{detail::lowerBounds(sliceData.logKFMin)};
    const std::array<double, 4> upperBounds{detail::upperBounds(sliceData.logKFMax)};

    if (!prevParams)
    {

        optimizer.setGuessBounds(coldGuess(), lowerBounds, upperBounds);

        return;
    }

    const Params<double> prevParamsD{prevParams->template as<double>()};
    optimizer.setGuessBounds(warmGuess(prevParamsD), lowerBounds, upperBounds);
}
} // namespace uv::models::svi::detail