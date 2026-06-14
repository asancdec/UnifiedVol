// SPDX-License-Identifier: Apache-2.0

#include "Math/Functions/Primitive.hpp"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <numbers>

TEST(MathPrimitive, NormalPdfAndCdfHaveKnownValues)
{
    EXPECT_NEAR(uv::math::normalPDF(0.0), 1.0 / std::sqrt(2.0 * std::numbers::pi), 1e-15);
    EXPECT_NEAR(uv::math::normalCDF(0.0), 0.5, 1e-15);
    EXPECT_NEAR(uv::math::normalCDF(1.0), 0.8413447460685429, 1e-15);
    EXPECT_NEAR(uv::math::normalCDF(-1.0), 0.15865525393145707, 1e-15);
}

TEST(MathPrimitive, ComplexHelpersMatchStandardLibrary)
{
    const std::complex<double> z{0.25, -0.2};

    const auto inv = uv::math::invComplex(z);
    const auto log1p = uv::math::log1pComplex(z);
    const auto expm1 = uv::math::expm1Complex(z);

    EXPECT_NEAR(std::real(inv), std::real(1.0 / z), 1e-15);
    EXPECT_NEAR(std::imag(inv), std::imag(1.0 / z), 1e-15);
    EXPECT_NEAR(std::real(log1p), std::real(std::log(1.0 + z)), 1e-15);
    EXPECT_NEAR(std::imag(log1p), std::imag(std::log(1.0 + z)), 1e-15);
    EXPECT_NEAR(std::real(expm1), std::real(std::exp(z) - 1.0), 1e-15);
    EXPECT_NEAR(std::imag(expm1), std::imag(std::exp(z) - 1.0), 1e-15);
}

TEST(MathPrimitive, Cosm1MatchesCosMinusOneNearZero)
{
    const double x = 1e-8;
    const double expected = -0.5 * x * x + (x * x * x * x) / 24.0;

    EXPECT_NEAR(uv::math::cosm1(x), expected, 1e-30);
}
