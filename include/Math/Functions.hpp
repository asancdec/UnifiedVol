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

#include "Utils/Types.hpp"
#include "Core/Matrix.hpp"

#include <concepts>

namespace uv::math
{
    /**
     * @brief Numerically stable complex log1p: log(1 + z).
     *
     * Implements the robust complex logarithm formulation from Appendix A of
     * Andersen–Lake, specialized to log(1 + z) for z = a + i b.
     *
     * The function uses a small-box branch when both real and imaginary parts
     * are small:
     *
     *   - If |a|, |b| < 1/2:
     *       Re(log(1 + z)) = 0.5 * log1p(a^2 + 2a + b^2)
     *       Im(log(1 + z)) = atan2(b, 1 + a)
     *
     *     This avoids catastrophic cancellation in log(1 + z) when z is close
     *     to zero.
     *
     *   - Otherwise:
     *       log(1 + z) = std::log(1 + z)
     *
     * @tparam T  Floating-point type (e.g. float, double).
     * @param z   Complex argument.
     *
     * @return log(1 + z) evaluated in a numerically robust way.
     *
     * ### REFERENCES
     * - L. Andersen, M. Lake (2018), *Robust High-Precision Option Pricing by
     *   Fourier Transforms: Contour Deformations and Double-Exponential
     *   Quadrature*, Bank of America Merrill Lynch, Appendix A.
     */
    template <std::floating_point T>
    Complex<T> log1pComplex(const Complex<T>& z) noexcept;

    /**
     * @brief Numerically stable evaluation of cos(b) - 1.
     *
     * Implements the cosm1 helper used in Appendix A of Andersen–Lake:
     *
     *   cosm1(b) = cos(b) - 1 = -2 sin²(b / 2),
     *
     * which replaces the naïve expression cos(b) - 1 that suffers severe
     * cancellation when |b| is small.
     *
     * @tparam T  Floating-point type.
     * @param b   Real argument.
     *
     * @return cos(b) - 1, computed in a cancellation-resistant form.
     *
     * ### REFERENCES
     * - L. Andersen, M. Lake (2018), *Robust High-Precision Option Pricing by
     *   Fourier Transforms: Contour Deformations and Double-Exponential
     *   Quadrature*, Appendix A.
     */
    template <std::floating_point T>
    T cosm1(T b) noexcept;

    /**
     * @brief Numerically stable complex expm1: exp(z) - 1.
     *
     * Implements the robust complex expm1 from Appendix A of Andersen–Lake:
     * for z = a + i b, use a small-|z| branch when |z| < 1 and fall back to
     * exp(z) - 1 otherwise.
     *
     * For |z| < 1:
     *
     *   - Use cosm1(b) := cos(b) - 1 (stable for small b),
     *   - Use em1 := expm1(a) = e^a - 1 (stable for small a),
     *
     *   and compute
     *
     *     Re(expm1(z)) = em1 * (cm1 + 1) + cm1
     *     Im(expm1(z)) = sin(b) * e^a
     *
     *   where cm1 = cosm1(b).
     *
     * For |z| >= 1:
     *
     *     expm1(z) = std::exp(z) - 1.
     *
     * This avoids catastrophic cancellation in exp(z) - 1 when z is close to
     * zero and is particularly important in Fourier-based option pricing
     * algorithms where many complex exponentials are evaluated near the origin.
     *
     * @tparam T  Floating-point type (e.g. float, double).
     * @param z   Complex argument.
     *
     * @return exp(z) - 1 evaluated in a numerically robust way.
     *
     * ### REFERENCES
     * - L. Andersen, M. Lake (2018), *Robust High-Precision Option Pricing by
     *   Fourier Transforms: Contour Deformations and Double-Exponential
     *   Quadrature*, Appendix A.
     */
    template <std::floating_point T>
    Complex<T> expm1Complex(const Complex<T>& z) noexcept;

    /**
     * @brief Standard normal cumulative distribution function Φ(x).
     *
     * Computes:
     *      Φ(x) = P(Z ≤ x),  Z ~ N(0, 1).
     *
     * @tparam T Floating-point type.
     * @param x  Evaluation point.
     * @return   Standard normal CDF at x.
     */
    template <std::floating_point T>
    T normalCDF(T x) noexcept;

    /**
     * @brief Standard normal probability density function φ(x).
     *
     * Computes:
     *      φ(x) = (1 / sqrt(2π)) * exp(-0.5 * x²).
     *
     * @tparam T Floating-point type.
     * @param x  Evaluation point.
     * @return   Standard normal PDF at x.
     */
    template <std::floating_point T>
    T normalPDF(T x) noexcept;

    /**
     * @brief Black–Scholes European option price.
     *
     * Computes the value of a European call or put option under
     * the Black–Scholes model with continuous dividend yield:
     *
     *   Call:
     *      C = S e^{-q t} Φ(d1) − K e^{-r t} Φ(d2)
     *
     *   Put:
     *      P = K e^{-r t} Φ(-d2) − S e^{-q t} Φ(-d1)
     *
     * @tparam T   Floating-point type.
     * @param t    Time to maturity (years).
     * @param r    Risk-free rate (continuous compounding).
     * @param q    Dividend yield (continuous compounding).
     * @param vol  Volatility σ.
     * @param S    Spot price.
     * @param K    Strike price.
     * @param isCall True for call, false for put.
     *
     * @return Option price (call or put).
     */
    template <std::floating_point T>
    T blackScholes(T t,
        T r,
        T q,
        T vol,
        T S,
        T K,
        bool isCall = true) noexcept;

     /**
     * @brief Black–Scholes pricing over a volatility surface.
     *
     * Prices European options across multiple maturities and strikes.
     * For each maturity index @p i, a full strike slice is priced using:
     *
     * @code
     *   prices[i][j] = BS(t[i], r[i], q[i], vol[i][j], S, K[i][j])
     * @endcode
     *
     * This overload internally calls the vector (slice) Black–Scholes
     * implementation for each maturity.
     *  
     * @param t      Vector of maturities (years).
     * @param r      Vector of risk-free rates.
     * @param q      Vector of dividend yields.
     * @param vol    Matrix of volatilities [maturity][strike].
     * @param S      Spot price.
     * @param K      Vector of strikes per maturity.
     * @param isCall True for calls, false for puts.
     *
     * @return Matrix of option prices [maturity][strike].
     *
     * @pre t.size() == r.size() == q.size() == vol.size()
     */
    core::Matrix<Real> blackScholes(const Vector<Real>& t,
        const Vector<Real>& r,
        const Vector<Real>& q,
        const core::Matrix<Real>& vol,
        Real S,
        const Vector<Real>& K,
        bool isCall = true);

    /**
     * @brief Black–Scholes Vega (∂Price / ∂σ).
     *
     * Computes:
     *      Vega = S e^{-q t} φ(d1) sqrt(t).
     *
     * Identical for calls and puts under the Black–Scholes model.
     *
     * @tparam T Floating-point type.
     * @param d1 Black–Scholes d1 term.
     * @param t  Time to maturity (years).
     * @param q  Dividend yield (continuous compounding).
     * @param S  Spot price.
     *
     * @return Vega in price units per unit volatility.
     */
    template <std::floating_point T>
    T vegaBS(T d1,
        T t,
        T q,
        T S) noexcept;

    /**
     * @brief Black–Scholes Volga (Vomma): ∂²Price / ∂σ².
     *
     * Computes:
     *      Volga = Vega * (d1 * d2 / σ),
     *   where d2 = d1 − σ sqrt(t).
     *
     * Identical for calls and puts under the Black–Scholes model.
     *
     * @tparam T   Floating-point type.
     * @param vega Precomputed Black–Scholes Vega.
     * @param d1   Black–Scholes d1 term.
     * @param t    Time to maturity (years).
     * @param vol  Volatility σ.
     *
     * @return Volga (second volatility derivative).
     */
    template <std::floating_point T>
    T volgaBS(T vega,
        T d1,
        T t,
        T vol) noexcept;

    /**
     * @brief Compute the Black–Scholes implied volatility using Halley's method.
     *
     * Finds the volatility such that the Black–Scholes price matches a given
     * call price. The function solves:
     *
     *   BS(vol) - callPrice = 0
     *
     * using Halley's iterative root-finding method, which combines Vega and
     * Volga for fast and stable convergence.
     *
     * @tparam T Floating-point type.
     *
     * @param callPrice  Market option price.
     * @param t          Time to maturity (years), must be positive.
     * @param r          Continuously compounded risk-free rate.
     * @param q          Continuous dividend yield.
     * @param S          Spot price of the underlying, must be positive.
     * @param K          Option strike, must be positive.
     *
     * @return Implied volatility such that the Black–Scholes price matches
     *         the given market price.
     *
     * @details
     * - The initial volatility guess is based on log-moneyness:
     *     vol0 = sqrt(2 * |log(F / K)| / t),
     *   where F = S * exp((r - q) * t).
     * - The algorithm performs up to 100 iterations with a tolerance of 1e-14
     *   on the pricing residual.
     *
     * @warning Emits a warning if convergence is not reached within the
     *          maximum number of iterations.
     */
    template <std::floating_point T>
    T impliedVolBS(T callPrice,
        T t,
        T r,
        T q,
        T S,
        T K);

    /**
     * @brief Compute Black–Scholes implied volatilities over a surface.
     *
     * Computes implied volatilities across multiple maturities using a common
     * strike grid. For each maturity index i, a full strike slice is solved via:
     *
     *   out[i] = impliedVolBS(callPrices[i], t[i], r[i], q[i], S, K, isCall)
     *
     * This overload internally calls the vector (slice) impliedVolBS implementation
     * for each maturity.
     *
     * @tparam T Floating-point type.
     *
     * @param callPrices Matrix of market option prices [maturity][strike].
     * @param t          Vector of maturities (years).
     * @param r          Vector of risk-free rates.
     * @param q          Vector of dividend yields.
     * @param S          Spot price.
     * @param K          Common strike vector used for all maturities.
     *
     * @return Matrix of implied volatilities [maturity][strike].
     *
     * @pre t.size() == r.size() == q.size() == callPrices.size()
     */
    template <std::floating_point T>
    core::Matrix<T> impliedVolBS(
        const core::Matrix<T>& callPrices,
        const Vector<T>& t,
        const Vector<T>& r,
        const Vector<T>& q,
        T S,
        const Vector<T>& K);

} // namespace uv::math

#include "Functions.inl"