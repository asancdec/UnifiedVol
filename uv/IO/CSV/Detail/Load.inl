// SPDX-License-Identifier: Apache-2.0

#include "Core/Generate.hpp"

namespace uv::io::csv::load
{

template <std::floating_point T> core::MarketState<T> marketState(
    const std::filesystem::path& path,
    const core::MarketData<T>& marketData,
    const detail::Options& opt
)
{
    auto [maturities, moneyness, vol] =
        detail::readLabeledMatrixCsv<T>(path.string(), opt);

    return ::uv::core::generateMarketState<T>(marketData, maturities, moneyness, vol);
}

} // namespace uv::io::csv::load
