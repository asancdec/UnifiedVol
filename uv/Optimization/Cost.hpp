// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <span>

namespace uv::opt::cost
{
template <std::floating_point T> struct WeightATM
{
    T wATM{1.0};
    T k0{1.0};
};

template <std::floating_point T> void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate = true
);
namespace detail
{

template <std::floating_point T> void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
);
}
} // namespace uv::opt::cost

#include "Optimization/Detail/Cost.inl"