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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under the License.
 */

#include <Models/SVI/Calibrate/Detail/Objective.hpp>

#include <cmath>

namespace uv::models::svi::detail
{

[[gnu::hot]] double
objectiveThunk(unsigned, const double* x, double* grad, void* data) noexcept
{
#if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
    const double* __restrict xx{x};
    double* __restrict gradOut{grad};
#else
    const double* xx{x};
    double* gradOut{grad};
#endif

    const auto& c{*static_cast<const ObjectiveContexts*>(data)};

    const double atm{c.atmTotalVariance};
    const std::size_t N{c.n};

    const double b{xx[0]};
    const double rho{xx[1]};
    const double mm{xx[2]};
    const double sigma{xx[3]};

    const double s2{sigma * sigma};

    const double R0{std::sqrt(mm * mm + s2)};
    const double invR0{1.0 / R0};

    const double brho{b * rho};
    const double bSigma{b * sigma};
    const double c0{atm - b * R0};

    const double* __restrict kptr{c.k};
    const double* __restrict wptr{c.wM};

    if (!gradOut)
    {
        double SSE{0.0};

        for (std::size_t i{0}; i < N; ++i)
        {
            const double k{*kptr++};
            const double wm{*wptr++};

            const double xi{k - mm};
            const double R{std::sqrt(xi * xi + s2)};

            const double wK{c0 + brho * k + b * R};

            const double r{wK - wm};
            SSE += r * r;
        }

        return SSE;
    }

    double SSE{0.0};
    double g0{0.0};
    double g1{0.0};
    double g2{0.0};
    double g3{0.0};

    const double mInvR0{mm * invR0};

    for (std::size_t i{0}; i < N; ++i)
    {
        const double k{*kptr++};
        const double wm{*wptr++};

        const double xi{k - mm};
        const double R{std::sqrt(xi * xi + s2)};
        const double invR{1.0 / R};

        const double wK{c0 + brho * k + b * R};

        const double r{wK - wm};
        SSE += r * r;

        const double twoR{2.0 * r};
        const double twoRB{twoR * b};

        g0 += twoR * (rho * k + R - R0);
        g1 += twoRB * k;
        g2 += -twoRB * (mInvR0 + xi * invR);
        g3 += twoR * bSigma * (invR - invR0);
    }

    gradOut[0] = g0;
    gradOut[1] = g1;
    gradOut[2] = g2;
    gradOut[3] = g3;

    return SSE;
}
} // namespace uv::models::svi::detail