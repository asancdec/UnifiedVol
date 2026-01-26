// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

#include "Core/Types.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <span>

namespace uv::core
{
/**
 * @brief Generate an evenly spaced 1D grid.
 *
 * Returns an array of values from bound1 to bound2 (inclusive).
 */
template <std::floating_point T, std::size_t N>
constexpr std::array<T, N> generateGrid(T bound1, T bound2) noexcept;

/**
 * @brief Evaluate a function on an array
 *
 * Returns a copy of the evaluated array
 */
template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept;

/**
 * @brief Evaluate a function in-place on an array
 *
 * Returns the same modified array
 */
template <std::floating_point T, std::size_t N, typename F>
void evalInplace(std::array<T, N>& grid, F&& f) noexcept;

/**
 * @brief Accumulate the elements in a span
 *
 * Returns the sum of the elements
 */
template <std::floating_point T> T sum(std::span<const T> x) noexcept;

/**
 * @brief Multiply all elements of a vector by a scalar.
 */
template <std::floating_point T> Vector<T> multiply(std::span<const T> v, T x) noexcept;

/**
 * @brief Compute element wise reciprocal of a vector.
 */
template <std::floating_point T> Vector<T> reciprocal(std::span<const T> v) noexcept;

/**
 * @brief Element wise multiplication of two vectors.
 */
template <std::floating_point T>
Vector<T> hadamard(std::span<const T> a, std::span<const T> b);

/**
 * @brief Generate a sequence of consecutive values.
 *
 * Returns a vector of length n starting from start.
 */
template <typename T> Vector<T> makeSequence(std::size_t n, T start) noexcept;

/**
 * @brief Return the minimum element of a vector.
 *
 * The input must be non empty.
 */
template <typename T> T minValue(std::span<const T> x);

/**
 * @brief Return the maximum element of a vector.
 *
 * The input must be non empty.
 */
template <typename T> T maxValue(std::span<const T> x);

/**
 * @brief Convert a vector to another value type.
 *
 * Returns a new vector with each element static casted.
 */
template <typename To, typename From>
Vector<To> convertVector(const Vector<From>& x) noexcept;

} // namespace uv::core

#include "Functions.inl"