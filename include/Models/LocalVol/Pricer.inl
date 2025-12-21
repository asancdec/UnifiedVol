// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.inl
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

#include "Math/PDE/Functions.hpp"
#include "Math/Integration/Functions.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <iomanip>
#include <utility>

namespace uv::models::localvol
{
    template
        <
        std::floating_point T,
        std::size_t nT,
        std::size_t nX,
        typename Fn
        >
    Pricer<T, nT, nX, Fn>::Pricer(T t,
        T r,
        T F,
        T K,
        Fn payoff,
        const std::array<T, nX>& xGrid)
        :
        t_(t),
        r_(r),
        F_(F),
        K_(K),
        payoff_(std::move(payoff)),
        xGrid_(xGrid),
        dt_(t_ / T(nT - 1)),
        dx_(xGrid_[1] - xGrid_[0]),
        pdfGrid_(math::pde::fokkerPlanckInit<T, nX>(T{ 0 }, xGrid_))
    {
        // TODO: Validate inputs
        // TODO: Set grid of strikes and tenors to price
    }

    template
        <
        std::floating_point T,
        std::size_t nT,
        std::size_t nX,
        typename Fn
        >
    T Pricer<T, nT, nX, Fn>::price(T localVar)
    {
        // TODO: Assume non constant localVar
        // TODO:: derivative of vol

        B_.fill(T{0.5} * localVar);
        C_.fill(T{ 0.5 } *localVar);

        // ---------- Solve Fokker-Plank PDE ----------

        math::pde::fokkerPlanckLog<Real, nT, nX>
        (
            pdfGrid_,
            B_,
            C_,
            dt_,
            dx_
        );

        // ---------- Price the payoff ----------

        return std::exp(-r_ * 3.0) *
                math::integration::trapezoidalWeighted
                (
                    [this] (T x) ->T
                    {
                        return payoff_(K_, normalizeForward_(F_) * std::exp(x));
                    },
                    pdfGrid_,
                    xGrid_
                );
        }

    template
        <
        std::floating_point T,
        std::size_t nT,
        std::size_t nX,
        typename Fn
        >
    T Pricer<T, nT, nX, Fn>::normalizeForward_(T F) const noexcept
    {
        return F / math::integration::trapezoidalWeighted
        (
            [](T x) {return std::exp(x); },
            pdfGrid_,
            xGrid_
        );
    }


} // namespace uv::models::localvol