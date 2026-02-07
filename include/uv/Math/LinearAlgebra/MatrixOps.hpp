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

#include "Core/Matrix.hpp"

#include <concepts>
#include <functional>

namespace uv::math::linear_algebra
{
template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t>
core::Matrix<T> generateIndexed(std::size_t rows, std::size_t cols, F&& f);

template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t, T>
core::Matrix<T> transformIndexed(const core::Matrix<T>& m, F&& f);

template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t, T>
void transformIndexedInplace(core::Matrix<T>& m, F&& f);

template <std::floating_point T>
core::Matrix<T> hadamard(const core::Matrix<T>& lhs, const core::Matrix<T>& rhs);

template <std::floating_point T>
void hadamardInplace(core::Matrix<T>& lhs, const core::Matrix<T>& rhs);

template <std::floating_point T>
core::Matrix<T> hadamard(const core::Matrix<T>& lhs, const Vector<T>& rhs);

template <std::floating_point T>
void hadamardInplace(core::Matrix<T>& lhs, const Vector<T>& rhs);

template <std::floating_point T>
core::Matrix<T> divide(const core::Matrix<T>& lhs, const core::Matrix<T>& rhs);

template <std::floating_point T>
void divideInplace(core::Matrix<T>& lhs, const core::Matrix<T>& rhs);

template <std::floating_point T>
core::Matrix<T> reciprocal(const core::Matrix<T>& m) noexcept;

template <std::floating_point T> void reciprocalInplace(core::Matrix<T>& m) noexcept;

template <std::floating_point T>
core::Matrix<T> square(const core::Matrix<T>& m) noexcept;

template <std::floating_point T> void squareInplace(core::Matrix<T>& m) noexcept;

template <std::floating_point T> core::Matrix<T> sqrt(const core::Matrix<T>& m);

template <std::floating_point T> void sqrtInplace(core::Matrix<T>& m);

template <typename To, typename From>
core::Matrix<To> convertMatrix(const core::Matrix<From>& A) noexcept;

} // namespace uv::math::linear_algebra

#include "Math/LinearAlgebra/Detail/MatrixOps.inl"