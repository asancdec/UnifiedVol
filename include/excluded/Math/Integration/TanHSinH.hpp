// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

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

    Node generateNode(T nh) const noexcept;

  public:
    TanHSinH();

    template <typename F> T integrateZeroToInf(F&& f) const noexcept;

    template <std::size_t M, typename F>
    std::array<T, M> integrateZeroToInfMulti(F&& f) const noexcept;

    void printGrid() const noexcept;
};
} // namespace uv::math::integration

#include "TanHSinH.inl"
