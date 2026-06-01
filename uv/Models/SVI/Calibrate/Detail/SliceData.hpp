// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <span>

namespace uv::models::svi::detail
{
struct SliceData
{
    const double atmTotalVariance;
    const double logKFMin;
    const double logKFMax;

    explicit SliceData(
        std::span<const double> logKF,
        std::span<const double> totalVariance
    );
};
} // namespace uv::models::svi::detail