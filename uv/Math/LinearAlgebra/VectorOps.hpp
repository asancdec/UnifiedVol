// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <span>

namespace uv::math::linear_algebra
{

template <std::floating_point T, std::size_t N, typename F>
std::array<T, N> eval(std::array<T, N> grid, F&& f) noexcept;

template <std::floating_point T, std::size_t N, typename F>
void evalInplace(std::array<T, N>& grid, F&& f) noexcept;

template <std::floating_point T> T sum(std::span<const T>) noexcept;

template <std::floating_point T> Vector<T> multiply(std::span<const T> v, T x);

template <std::floating_point T> Vector<T> reciprocal(std::span<const T>);

template <std::floating_point T> Vector<T> exponential(std::span<const T>);

template <std::floating_point T> Vector<T> squareRoot(std::span<const T>);

template <std::floating_point T> void squareRootInplace(std::span<T>, std::span<const T>);

template <std::floating_point T>
Vector<T> hadamard(std::span<const T> a, std::span<const T> b);

template <typename T> Vector<T> makeSequence(std::size_t n, T start);

template <typename T> T minValue(std::span<const T>);

template <typename T> T maxValue(std::span<const T>);

} // namespace uv::math::linear_algebra

#include "Math/LinearAlgebra/Detail/VectorOps.inl"
