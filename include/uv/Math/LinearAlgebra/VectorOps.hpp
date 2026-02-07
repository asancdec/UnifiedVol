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

#include "Base/Alias.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <span>

namespace uv::math::linear_algebra
{

template <std::floating_point T, std::size_t N>
constexpr std::array<T, N> generateGrid(T bound1, T bound2) noexcept;

template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept;

template <std::floating_point T, std::size_t N, typename F>
void evalInplace(std::array<T, N>& grid, F&& f) noexcept;

template <std::floating_point T> T sum(std::span<const T> x) noexcept;

template <std::floating_point T> Vector<T> multiply(std::span<const T> v, T x) noexcept;

template <std::floating_point T> Vector<T> reciprocal(std::span<const T> v) noexcept;

template <std::floating_point T>
Vector<T> hadamard(std::span<const T> a, std::span<const T> b);

template <typename T> Vector<T> makeSequence(std::size_t n, T start) noexcept;

template <typename T> T minValue(std::span<const T> x);

template <typename T> T maxValue(std::span<const T> x);

template <typename To, typename From>
Vector<To> convertVector(const Vector<From>& x) noexcept;

} // namespace uv::math::linear_algebra

#include "Math/LinearAlgebra/Detail/VectorOps.inl"