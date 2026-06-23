// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Models/SVI/Params.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include <concepts>
#include <span>

namespace uv::models::svi
{
template <std::floating_point T, opt::nlopt::Algorithm Algo> Vector<Params<T>> calibrate(
    const core::VolSurface<T>& volSurface,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    bool printParams = false
);

template <std::floating_point T, opt::nlopt::Algorithm Algo> Vector<Params<T>> calibrate(
    std::span<const T> maturities,
    const core::Matrix<T>& logKF,
    const core::Matrix<T>& totalVariance,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    bool printParams = false
);
} // namespace uv::models::svi

#include "Models/SVI/Calibrate/Detail/Calibrate.inl"
