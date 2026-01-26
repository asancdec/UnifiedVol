// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Public interface and documentation for SVI calibration and surface
 * construction.
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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "Core/Matrix/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Optimization/NLopt/Optimizer.hpp"
#include "Models/SVI/Params.hpp"

#include <array>
#include <concepts>
#include <nlopt.hpp>
#include <span>
#include <vector>

namespace opt = uv::math::opt;

namespace uv::models::svi
{

/**
 * @brief Calibrate an SVI surface from a volatility surface.
 *
 * Extracts maturities and total variances from a @ref core::VolSurface
 * and performs slice-by-slice SVI calibration.
 *
 * @note
 * Calibration follows the arbitrage-free methodology of
 * Tahar Ferhati (2020), *Robust Calibration for SVI Model Arbitrage Free*,
 * enforcing butterfly and calendar-spread constraints.
 *
 * @tparam T     Floating-point type.
 * @tparam Algo  NLopt algorithm type.
 *
 * @param volSurface Market volatility surface.
 * @param prototype  Prototype NLopt optimizer.
 *
 * @return Calibrated SVI parameters, one per maturity.
 */
template <std::floating_point T, ::nlopt::algorithm Algo>
Vector<Params<T>> calibrate(
    const core::VolSurface<T>& volSurface,
    const opt::nlopt::Optimizer<4, Algo>& prototype
);

/**
 * @brief Calibrate SVI parameters independently per maturity.
 *
 * Fits one SVI slice per tenor using total variance data.
 *
 * @note
 * Arbitrage-free constraints follow Ferhati (2020),
 * *Robust Calibration for SVI Model Arbitrage Free*.
 *
 * @tparam T     Floating-point type.
 * @tparam Algo  NLopt algorithm type.
 *
 * @param tenors    Maturities (years).
 * @param kMatrix   Log-forward moneyness grid.
 * @param wMMatrix  Market total variance matrix.
 * @param prototype Prototype NLopt optimizer.
 *
 * @return Calibrated SVI parameters, one per maturity.
 */
template <std::floating_point T, ::nlopt::algorithm Algo>
Vector<Params<T>> calibrate(
    const Vector<T>& tenors,
    const core::Matrix<T>& kMatrix,
    const core::Matrix<T>& wMMatrix,
    const opt::nlopt::Optimizer<4, Algo>& prototype
);

/**
 * @brief Build an SVI-implied total variance surface.
 *
 * Uses calibrated SVI parameters to compute total variance
 * on the same grid as an existing volatility surface.
 *
 * The returned surface shares the same strikes, tenors,
 * and market metadata as the input surface.
 *
 * @param volSurface Reference volatility surface defining the grid.
 * @param params     Calibrated SVI parameters (one per tenor).
 *
 * @return New volatility surface with SVI total variance.
 *
 * @throws ErrorCode::InvalidArgument if parameter dimensions mismatch.
 */
template <std::floating_point T>
core::VolSurface<T>
buildSurface(const core::VolSurface<T>& volSurface, const Vector<Params<T>>& params);

/**
 * @brief SVI convexity function g(k) used for butterfly arbitrage checks.
 *
 * Computes the SVI convexity condition function g(k) for a single slice.
 * In no-butterfly-arbitrage SVI calibration, the constraint g(k) >= 0 is
 * enforced on a strike grid.
 *
 * @param a     SVI level parameter.
 * @param b     SVI slope/amplitude parameter.
 * @param rho   SVI skew parameter (-1 < rho < 1).
 * @param m     SVI horizontal shift parameter.
 * @param sigma SVI curvature/smoothness parameter (sigma > 0).
 * @param k     Log-forward moneyness.
 *
 * @return Value of g(k).
 */
template <std::floating_point T>
[[gnu::hot]] T gk(T a, T b, T rho, T m, T sigma, T k) noexcept;

namespace detail
{
/// @brief Calendar-spread constraint context (internal).
struct CalendarContexts;

/// @brief Objective function context (internal).
struct ObjectiveContexts;

/// @brief Precomputed values for g(k) and its gradient (internal).
struct GkCache;

/// @brief Convexity constraint context for a single k (internal).
struct ConvexityContexts;

/**
 * @brief Validate calibration input grids and dimensions.
 *
 * Checks that:
 * - tenors is non-empty and strictly increasing
 * - kMatrix and wMMatrix have the same number of rows as tenors
 * - matrices are rectangular
 * - each kSlice is strictly increasing
 * - kMatrix and wMMatrix dimensions match
 *
 * @param tenors  Maturities (years).
 * @param kMatrix Log-forward moneyness grid per tenor.
 * @param wMMatrix Market total variance per tenor and strike.
 *
 * @throws ErrorCode::InvalidArgument if validation fails.
 */
template <std::floating_point T>
void validateInputs(
    const Vector<T>& tenors,
    const core::Matrix<T>& kMatrix,
    const core::Matrix<T>& wMMatrix
);

/**
 * @brief Calibrate one SVI slice for a single maturity.
 *
 * Runs an NLopt optimization for one tenor using:
 * - positive total variance constraint
 * - Roger–Lee slope constraints
 * - convexity constraints g(k) >= 0 on the provided k grid
 * - optional calendar constraint versus the previous slice
 *
 * The ATM total variance is used to pin the `a` parameter via aParam().
 *
 * @tparam Algo NLopt algorithm.
 *
 * @param t          Slice maturity (years).
 * @param kSlice     Log-forward moneyness grid for this tenor (strictly
 * increasing).
 * @param wKSlice    Market total variance values at kSlice.
 * @param prototype  Prototype optimizer used to create a fresh instance.
 * @param prevParams Previous slice parameters (nullptr for first slice).
 * @param numStrikes Number of strikes in the slice.
 *
 * @return Calibrated SVI parameters for this slice.
 */
template <std::floating_point T, ::nlopt::algorithm Algo>
Params<T> calibrateSlice(
    T t,
    std::span<const double> kSlice,
    std::span<const double> wKSlice,
    const opt::nlopt::Optimizer<4, Algo>& prototype,
    const Params<T>* prevParams,
    std::size_t numStrikes
);

/**
 * @brief Cold-start initial guess for (b, rho, m, sigma).
 *
 * Used when no previous calibration is available.
 */
std::array<double, 4> coldGuess() noexcept;

/**
 * @brief Warm-start initial guess from previous parameters.
 *
 * Extracts (b, rho, m, sigma) from an existing parameter set.
 *
 * @param params Previously calibrated parameters.
 */
template <std::floating_point T>
std::array<double, 4> warmGuess(const Params<T>& params) noexcept;

/**
 * @brief Lower bounds for (b, rho, m, sigma).
 *
 * Bounds may depend on the minimum k in the slice to avoid pathological
 * parameter combinations and improve numerical stability.
 *
 * @param logKFMin Minimum log-forward moneyness in the slice.
 *
 * @return Lower bounds vector [b, rho, m, sigma].
 */
std::array<double, 4> lowerBounds(double logKFMin) noexcept;

/**
 * @brief Upper bounds for (b, rho, m, sigma).
 *
 * Bounds may depend on the maximum k in the slice to avoid pathological
 * parameter combinations and improve numerical stability.
 *
 * @param logKFMax Maximum log-forward moneyness in the slice.
 *
 * @return Upper bounds vector [b, rho, m, sigma].
 */
std::array<double, 4> upperBounds(double logKFMax) noexcept;

/**
 * @brief Add minimum total variance constraint.
 *
 * Enforces w(k) >= 0 (typically at a key location such as k=0 or across
 * a small set of points, depending on implementation).
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 */
template <::nlopt::algorithm Algo>
void addWMinConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept;

/**
 * @brief Add Roger–Lee minimum wing slope constraints.
 *
 * Enforces lower bounds on the asymptotic wing slopes implied by SVI
 * parameters to satisfy no-arbitrage requirements.
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 */
template <::nlopt::algorithm Algo>
void addMinSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept;

/**
 * @brief Add Roger–Lee maximum wing slope constraints.
 *
 * Enforces upper bounds on the asymptotic wing slopes implied by SVI
 * parameters to satisfy no-arbitrage requirements.
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 */
template <::nlopt::algorithm Algo>
void addMaxSlopeConstraint(opt::nlopt::Optimizer<4, Algo>& optimizer) noexcept;

/**
 * @brief Build a calendar-spread constraint context at strike k.
 *
 * Computes w_prev(k) from the supplied previous-slice SVI parameters and
 * returns a CalendarContexts entry used by calendar monotonicity constraints.
 *
 * @tparam T Floating-point type of the SVI parameters.
 * @param k Log-moneyness at which the constraint is enforced.
 * @param params SVI parameters of the previous maturity slice.
 * @param eps Non-negative slack for the calendar constraint.
 * @param atmWK ATM total variance (used for normalization/weighting).
 * @return Calendar constraint context for strike k.
 */
template <std::floating_point T>
CalendarContexts
genCalendarContext(double k, const Params<T>& params, double eps, double atmWK) noexcept;

/**
 * @brief Populate calendar-spread constraint contexts for a maturity slice.
 *
 * Fills @p out with one CalendarContexts entry per strike in @p kSlice and
 * optionally adds two extra points outside the market range using @p delta.
 *
 * @tparam T Floating-point type of the SVI parameters.
 * @tparam Algo NLopt algorithm.
 * @param out Output container for calendar constraint contexts (must outlive
 * optimize()).
 * @param numStrikes Number of strikes in @p kSlice.
 * @param params SVI parameters of the previous maturity slice.
 * @param optimizer Optimizer instance (for eps/slack settings).
 * @param kSlice Strike grid (log-moneyness), assumed sorted.
 * @param atmWK ATM total variance (used for normalization/weighting).
 * @param logKFMin Minimum log(K/F) of the market range.
 * @param logKFMax Maximum log(K/F) of the market range.
 * @param delta Extension added to both sides of the market range.
 */
template <std::floating_point T, ::nlopt::algorithm Algo>
void fillCalendarContexts(
    Vector<CalendarContexts>& calendarContexts,
    std::size_t numStrikes,
    const Params<T>& params,
    const opt::nlopt::Optimizer<4, Algo>& optimizer,
    std::span<const double> kSlice,
    double atmWK,
    double logKFMin,
    double logKFMax,
    double delta = 0.15
) noexcept;

/**
 * @brief Add calendar-spread constraint against a previous slice.
 *
 * Enforces total variance monotonicity in maturity:
 * w_current(k) >= w_previous(k) for each k in the grid.
 *
 * @note The contexts in `ctx` must remain alive until the optimization ends.
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 * @param ctx Vector of per-strike calendar contexts (one per k).
 */
template <::nlopt::algorithm Algo>
void addCalendarConstraint(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    Vector<CalendarContexts>& ctx
) noexcept;

/**
 * @brief Add convexity (no-butterfly) constraint at one k.
 *
 * Adds a single inequality constraint enforcing g(k) >= 0 at a specific
 * log-forward moneyness value.
 *
 * @note The context `ctx` must remain alive until the optimization ends.
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 * @param ctx Convexity context for one k.
 */
template <::nlopt::algorithm Algo>
void addConvexityConstraint(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    ConvexityContexts& ctx
) noexcept;

/**
 * @brief Set the least-squares objective for one slice.
 *
 * Sets the objective function to fit SVI total variance to the market
 * total variance on the provided k grid. An analytical gradient is used.
 *
 * @note The context `ctx` must remain alive until the optimization ends.
 *
 * @tparam Algo NLopt algorithm.
 * @param optimizer Optimizer instance.
 * @param ctx Objective context (market data pointers, sizes, ATM total
 * variance).
 */
template <::nlopt::algorithm Algo>
void setMinObjective(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const ObjectiveContexts& ctx
) noexcept;

/**
 * @brief Compute SVI total variance w(k).
 *
 * Evaluates the SVI total variance function at log-forward moneyness k.
 *
 * @param a     SVI level parameter.
 * @param b     SVI slope/amplitude parameter.
 * @param rho   SVI skew parameter.
 * @param m     SVI horizontal shift parameter.
 * @param sigma SVI curvature/smoothness parameter.
 * @param k     Log-forward moneyness.
 *
 * @return Total variance w(k).
 */
template <std::floating_point T>
[[gnu::hot]] T calculateWk(T a, T b, T rho, T m, T sigma, T k) noexcept;

/**
 * @brief Compute the pinned `a` parameter from ATM total variance.
 *
 * Computes `a` such that w(0) equals the provided ATM total variance.
 * Used when `a` is not optimized directly.
 *
 * @param atmWK ATM total variance w(0).
 * @param b     SVI slope/amplitude parameter.
 * @param rho   SVI skew parameter.
 * @param m     SVI horizontal shift parameter.
 * @param sigma SVI curvature/smoothness parameter.
 *
 * @return The implied `a` value for this slice.
 */
[[gnu::hot]] double
aParam(double atmWK, double b, double rho, double m, double sigma) noexcept;

/**
 * @brief Compute g(k) using precomputed cache values (internal).
 *
 * Fast evaluation of g(k) for convexity constraints using values stored
 * in @ref GkCache.
 *
 * @param p Precomputed cache for the current parameters and k.
 *
 * @return Value of g(k).
 */
[[gnu::hot]] double gk(const GkCache& p) noexcept;

/**
 * @brief Gradient of g(k) with respect to (b, rho, m, sigma).
 *
 * Computes the analytical gradient used by NLopt for convexity constraints.
 * The `a` parameter is assumed to be ATM-pinned and is handled consistently
 * via the cache / chain rule.
 *
 * @param b     SVI slope/amplitude parameter.
 * @param rho   SVI skew parameter.
 * @param m     SVI horizontal shift parameter.
 * @param sigma SVI curvature/smoothness parameter.
 * @param k     Log-forward moneyness.
 * @param p     Precomputed cache for the current parameters and k.
 *
 * @return Gradient vector [db, drho, dm, dsigma].
 */
[[gnu::hot]] std::array<double, 4>
gkGrad(double b, double rho, double m, double sigma, double k, const GkCache& p) noexcept;

} // namespace detail
} // namespace uv::models::svi

#include "Functions.inl"
