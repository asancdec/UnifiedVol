// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace uv::tests::tolerance
{
inline constexpr double Golden = 1e-8;
inline constexpr double FiniteDifference = 2e-4;
inline constexpr double InterpolationInvariant = 1e-14;
inline constexpr double PricingInvariant = 1e-13;
inline constexpr double NoArb = 5e-10;
inline constexpr double RoundTrip = 1e-10;
} // namespace uv::tests::tolerance
