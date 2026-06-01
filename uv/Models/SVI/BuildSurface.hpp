// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/Calibrate/Config.hpp"
#include "Models/SVI/Params.hpp"

#include <concepts>

namespace uv::models::svi
{
template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::MarketState<T>& marketState, const Config& config = {});

template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Config& config = {});

template <std::floating_point T> core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Vector<Params<T>>& params);

} // namespace uv::models::svi

#include "Models/SVI/Detail/BuildSurface.inl"
