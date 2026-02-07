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

#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"

#include "concepts"
#include <cstddef>
#include <span>

namespace uv::core
{

template <std::floating_point T> class VolSurface
{
  private:
    Vector<T> maturities_;
    std::size_t numMaturities_;

    Vector<T> strikes_;
    std::size_t numStrikes_;

    Vector<T> forwards_;
    Vector<T> moneyness_;

    Matrix<T> vol_;

  public:
    VolSurface() = delete;

    explicit VolSurface(
        std::span<const T> maturities,
        std::span<const T> forwards,
        std::span<const T> strikes,
        std::span<const T> moneyenss,
        const Matrix<T>& vol
    );

    std::size_t numMaturities() const noexcept;
    std::size_t numStrikes() const noexcept;

    std::span<const T> maturities() const noexcept;
    std::span<const T> forwards() const noexcept;
    std::span<const T> strikes() const noexcept;
    std::span<const T> moneyness() const noexcept;

    const Matrix<T>& vol() const noexcept;
};
} // namespace uv::core

#include "Core/Detail/VolSurface.inl"
