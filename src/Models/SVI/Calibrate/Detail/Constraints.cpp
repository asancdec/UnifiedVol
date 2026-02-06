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

#include <Models/SVI/Calibrate/Detail/Constraints.hpp>

#include <cmath>

namespace uv::models::svi::detail
{

[[gnu::hot]]
void calendarMConstraint(
    unsigned m,
    double* result,
    unsigned,
    const double* x,
    double* grad,
    void* data
) noexcept
{
#if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
    auto* __restrict res{result};
    auto* __restrict gradOut{grad};
    const double* __restrict xx{x};
#else
    auto* res{result};
    auto* gradOut{grad};
    const double* xx{x};
#endif

    const auto& ctx{*static_cast<const CalendarMContext*>(data)};

    constexpr unsigned nParams{4};

    const double atm{ctx.atmTotalVariance};
    const double eps{ctx.eps};

    const double b{xx[0]};
    const double rho{xx[1]};
    const double mm{xx[2]};
    const double sigma{xx[3]};

    const double s2{sigma * sigma};
    const double R0{std::sqrt(mm * mm + s2)};

    const double* __restrict kptr{ctx.logKF.data()};
    const double* __restrict pptr{ctx.prevWk.data()};

    if (!gradOut)
    {
        const double brho{b * rho};
        const double c{eps - atm + b * R0};

        for (unsigned i{0}; i < m; ++i)
        {
            const double k{*kptr++};
            const double prev{*pptr++};

            const double xi{k - mm};
            const double Rk{std::sqrt(xi * xi + s2)};

            *res++ = prev + c - b * Rk - brho * k;
        }
        return;
    }

    const double invR0{1.0 / R0};
    const double mInvR0{mm * invR0};
    const double bSigma{b * sigma};
    const double c0{eps - atm};

    double* __restrict gptr{gradOut};

    for (unsigned i{0}; i < m; ++i)
    {
        const double k{*kptr++};
        const double prev{*pptr++};

        const double xi{k - mm};
        const double Rk{std::sqrt(xi * xi + s2)};
        const double invRk{1.0 / Rk};

        const double dR{R0 - Rk};
        const double rhok{rho * k};

        *res++ = prev + c0 + b * (dR - rhok);

        gptr[0] = dR - rhok;
        gptr[1] = -(b * k);
        gptr[2] = b * (mInvR0 + xi * invRk);
        gptr[3] = bSigma * (invR0 - invRk);

        gptr += nParams;
    }
}

[[gnu::hot]]
void convexityMConstraint(
    unsigned m,
    double* result,
    unsigned,
    const double* x,
    double* grad,
    void* data
) noexcept
{
#if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
    auto* __restrict res{result};
    auto* __restrict gradOut{grad};
    const double* __restrict xx{x};
#else
    auto* res{result};
    auto* gradOut{grad};
    const double* xx{x};
#endif

    const auto& ctx{*static_cast<const ConvexityMContext*>(data)};

    constexpr double half{0.5};
    constexpr double quarter{0.25};
    constexpr unsigned nParams{4};

    const double atm{ctx.atmTotalVariance};

    const double b{xx[0]};
    const double rho{xx[1]};
    const double mm{xx[2]};
    const double sigma{xx[3]};

    const double s2{sigma * sigma};
    const double mm2{mm * mm};

    const double R0{std::sqrt(mm2 + s2)};
    const double invR0{1.0 / R0};

    const double a{atm + (b * rho) * mm - b * R0};

    const double rhoMmMinusR0{rho * mm - R0};
    const double rhoMinusMmInvR0{rho - mm * invR0};
    const double bSigma{b * sigma};

    const double* __restrict kptr{ctx.logKF.data()};

    if (!gradOut)
    {
        for (unsigned i{0}; i < m; ++i)
        {
            const double k{*kptr++};
            const double xkm{k - mm};

            const double R{std::sqrt(xkm * xkm + s2)};
            const double invR{1.0 / R};

            const double invR2{invR * invR};
            const double invR3{invR2 * invR};

            const double s2InvR3{s2 * invR3};

            const double t{rho + xkm * invR};
            const double rhoXkmR{rho * xkm + R};

            const double w{a + b * rhoXkmR};
            const double wInv{1.0 / w};

            const double w1{b * t};
            const double w1Sq{w1 * w1};

            const double A{1.0 - half * k * (w1 * wInv)};
            const double B{wInv + quarter};

            const double gval{(A * A) - quarter * w1Sq * B + half * (b * s2InvR3)};

            res[i] = -gval;
        }
        return;
    }

    double* __restrict gptr{gradOut};

    for (unsigned i{0}; i < m; ++i)
    {
        const double k{*kptr++};
        const double xkm{k - mm};

        const double R{std::sqrt(xkm * xkm + s2)};
        const double invR{1.0 / R};

        const double invR2{invR * invR};
        const double invR3{invR2 * invR};

        const double s2InvR3{s2 * invR3};
        const double s2InvR5{s2InvR3 * invR2};

        const double t{rho + xkm * invR};
        const double rhoXkmR{rho * xkm + R};

        const double w{a + b * rhoXkmR};
        const double wInv{1.0 / w};
        const double wInv2{wInv * wInv};

        const double w1{b * t};
        const double w1Sq{w1 * w1};

        const double A{1.0 - half * k * (w1 * wInv)};
        const double B{wInv + quarter};

        const double gval{(A * A) - quarter * w1Sq * B + half * (b * s2InvR3)};

        res[i] = -gval;

        const double Ak{A * k};

        const double dgdw{wInv2 * w1 * (Ak + quarter * w1)};

        const double dgdw1{(-Ak) * wInv - half * w1 * B};

        const double dg0{dgdw * (rhoXkmR + rhoMmMinusR0) + dgdw1 * t + half * s2InvR3};

        const double dg1{b * (dgdw * k + dgdw1)};

        const double dg2{
            b * (dgdw * (-t + rhoMinusMmInvR0) - dgdw1 * s2InvR3 + 1.5 * xkm * s2InvR5)
        };

        const double dg3{
            bSigma *
            (dgdw * (invR - invR0) - dgdw1 * (xkm * invR3) + (invR3 - 1.5 * s2InvR5))
        };

        gptr[0] = -dg0;
        gptr[1] = -dg1;
        gptr[2] = -dg2;
        gptr[3] = -dg3;
        gptr += nParams;
    }
}
} // namespace uv::models::svi::detail