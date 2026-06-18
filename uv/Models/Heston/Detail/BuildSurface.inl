// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Core/Generate.hpp"
#include "Core/Matrix.hpp"
#include "Math/Functions/Volatility.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/Heston/Calibrate/Config.hpp"

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N> core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const calibrate::Config& config
)
{
    price::Pricer<T, N> pricer{};

    pricer.setParams(calibrate::calibrate(volSurface, curve, config, pricer));
    return buildSurface(volSurface, curve, pricer);
}

template <std::floating_point T, std::size_t N> core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::MarketState<T>& marketState,
    const calibrate::Config& config
)
{
    return buildSurface<T, N>(volSurface, marketState.interestCurve, config);
}

template <std::floating_point T, std::size_t N> core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const price::Pricer<T, N>& pricer
)
{
    return core::generateVolSurface(
        volSurface,
        math::vol::impliedVol(pricer.callPrice(volSurface, curve), volSurface, curve)
    );
}

template <std::floating_point T, std::size_t N> core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::MarketState<T>& marketState,
    const price::Pricer<T, N>& pricer
)
{
    return buildSurface(volSurface, marketState.interestCurve, pricer);
}
} // namespace uv::models::heston
