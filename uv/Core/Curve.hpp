// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>
#include <cstddef>
#include <span>

namespace uv::core
{

template <std::floating_point T> class Curve
{
  private:
    std::size_t numMaturities_;
    Vector<T> maturities_;
    Vector<T> discountFactors_;

  public:
    Curve() = delete;

    explicit Curve(T continuouslyCompoundedRate, std::span<const T> maturities);

    T interpolateDF(T maturity, bool doValidate = true) const;

    Vector<T> interpolateDF(std::span<const T> maturities, bool doValidate = true) const;
};
} // namespace uv::core

#include "Core/Detail/Curve.inl"