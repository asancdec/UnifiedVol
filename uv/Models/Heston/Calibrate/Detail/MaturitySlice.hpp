// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/Matrix.hpp"
#include "Optimization/Cost.hpp"

#include <cstddef>
#include <span>

namespace uv::models::heston::calibrate::detail
{

struct MaturitySlice
{
    double t;
    double dF;
    double F;
    Vector<double> K;
    Vector<double> vol;
    Vector<double> w;

    MaturitySlice() = delete;

    explicit MaturitySlice(std::size_t capacity);
};

Vector<MaturitySlice> makeSlices(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& vol,
    const opt::cost::WeightATM<double>& weightATM
);

void validateInputs(
    const std::span<const double> maturities,
    const std::span<const double> discountFactors,
    const std::span<const double> forwards,
    const std::span<const double> strikes,
    const core::Matrix<double>& vol
);

} // namespace uv::models::heston::calibrate::detail
