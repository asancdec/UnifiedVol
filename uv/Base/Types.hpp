// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <complex>
#include <concepts>
#include <span>
#include <vector>

namespace uv
{
template <typename T> using Vector = std::vector<T>;

template <std::floating_point T> using Complex = std::complex<T>;

template <typename To, typename From> Vector<To> convertVector(const Vector<From>&);

template <typename To, typename From> Vector<To> convertVector(std::span<const From>);
} // namespace uv

#include "Base/Detail/Types.inl"
