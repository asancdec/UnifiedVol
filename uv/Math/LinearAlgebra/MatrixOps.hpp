// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Core/Matrix.hpp"

#include <concepts>

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

template <std::floating_point T> core::Matrix<T> reciprocal(const core::Matrix<T>&);

template <std::floating_point T> void reciprocalInplace(core::Matrix<T>&) noexcept;

template <std::floating_point T> core::Matrix<T> square(const core::Matrix<T>&);

template <std::floating_point T> void squareInplace(core::Matrix<T>&) noexcept;

template <std::floating_point T> core::Matrix<T> sqrt(const core::Matrix<T>&);

template <std::floating_point T> void sqrtInplace(core::Matrix<T>&);

} // namespace uv::math::linear_algebra

#include "Math/LinearAlgebra/Detail/MatrixOps.inl"
