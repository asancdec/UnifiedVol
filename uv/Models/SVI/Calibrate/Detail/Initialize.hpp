// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/SVI/Calibrate/Detail/SliceData.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <array>
#include <concepts>

namespace uv::models::svi::detail
{

template <std::floating_point T, opt::nlopt::Algorithm Algo> void setGuessBounds(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const Params<T>* prevParams,
    const SliceData& sliceData
) noexcept;

std::array<double, 4> coldGuess() noexcept;

std::array<double, 4> warmGuess(const Params<double>& params) noexcept;

std::array<double, 4> lowerBounds(double logKFMin) noexcept;

std::array<double, 4> upperBounds(double logKFMax) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Initialize.inl"