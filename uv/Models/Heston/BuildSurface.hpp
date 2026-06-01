// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Curve.hpp"
#include "Core/VolSurface.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Price/Pricer.hpp"

#include <concepts>
#include <cstddef>

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N = calibrate::defaultNodes>
core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const calibrate::Config& config = {}
);

template <std::floating_point T, std::size_t N> core::VolSurface<T> buildSurface(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const price::Pricer<T, N>& pricer
);

} // namespace uv::models::heston

#include "Models/Heston/Detail/BuildSurface.inl"
