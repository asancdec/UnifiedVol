// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <concepts>
#include <cstddef>

namespace uv::math::integration
{

template <std::floating_point T, std::size_t N> class TanHSinH
{
  private:
    struct Node
    {
        T w;
        T y;
        T x;
        T factorRight;
        T inputRight;
        T factorLeft;
        T inputLeft;
    };

    const T h_;
    std::array<Node, N> nodes_;

    Node generateNode(T) const noexcept;

  public:
    TanHSinH();

    template <typename F> T integrateZeroToInf(F&&) const noexcept;

    template <std::size_t M, typename F>
    std::array<T, M> integrateZeroToInfMulti(F&&) const noexcept;
};
} // namespace uv::math::integration

#include "Math/Integration/Detail/TanHSinH.inl"
