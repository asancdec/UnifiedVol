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

#include <Base/Alias.hpp>

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

    T discountFactor(T maturity, bool doValidate = true) const;

    Vector<T> discountFactor(std::span<const T> maturities, bool doValidate = true) const;
};
} // namespace uv::core

#include <Core/Detail/Curve.inl>